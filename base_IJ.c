#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"

#define MIN(a, b) (a) < (b) ? (a) : (b)
#ifndef TEST_CNT
#define TEST_CNT 1
#endif

HYPRE_BigInt read_binary(void * buf, const char * filename, int size_of_elems, long long start, HYPRE_BigInt nums) {
    assert(size_of_elems == 4 || size_of_elems == 8);
    int my_pid; MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
    if (my_pid == 0) printf("reading binary from %s\n", filename);

    FILE * fp = fopen(filename, "rb");
    if (fp == NULL) {
        printf("cannot open %s \n", filename);
        return -1;
    }
    if (fseek(fp, size_of_elems * start, SEEK_SET) != 0) {
        printf("Error! cannot move file pointer to %d-th bytes\n", size_of_elems * start);
        fclose(fp);
        return -1;
    }
    HYPRE_BigInt ret;
    
    ret = fread(buf, size_of_elems, nums, fp);
    return ret;
}

void check_input(const HYPRE_Real * input, HYPRE_BigInt num) {
    for (HYPRE_BigInt i = 0; i < num; i++)
        // assert(input[i] == input[i]);
        assert(abs(input[i]) < 1e20);
}

int main(int argc, char * argv[])
{
    int my_pid, num_procs;
    HYPRE_Int glb_nblks = -1, ndof = -1;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);

    int cnt = 1;
    const char * case_name = argv[cnt++];
    HYPRE_Int gx = atoi(argv[cnt++]),
        gy = atoi(argv[cnt++]),
        gz = atoi(argv[cnt++]);
    if (strcmp(case_name, "LASER") == 0) {
        ndof = 1;
        glb_nblks = gx * gy * gz;
    }
    else if (strcmp(case_name, "LASER-3T") == 0) {
        ndof = 3;
        glb_nblks = gx * gy * gz;
    }
    else if (strcmp(case_name, "OIL-IMPEC") == 0) {
        ndof = 1;
        glb_nblks = gx * gy * gz + 5;
    }
    else if (strcmp(case_name, "GRAPES") == 0) {
        ndof = 1;
        glb_nblks = gx * gy * gz;
    }
    else if (strcmp(case_name, "OIL-FIM") == 0) {
        HYPRE_Int num_wells = 5;
        if (my_pid == 0) printf("OIL-FIM num wells: %d\n", num_wells);
        ndof = 4;
        glb_nblks = gx * gy * gz + num_wells;
    }
    else if (strcmp(case_name, "SOLID") == 0) {
        ndof = 3;
        glb_nblks = gx * gy * gz;
    }
    else {
        if (my_pid == 0) printf("Error: Unrecognized case name %s\n", case_name);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    HYPRE_Int loc_nrows, glb_nrows = glb_nblks * ndof;
    HYPRE_Int loc_nblks = -1;
    HYPRE_Int ilower, iupper;// 本进程负责的起始行号 和末尾行号（闭区间）
    
    if (my_pid == 0) printf("using blk-partition!\n");
    loc_nblks = glb_nblks / num_procs;
    HYPRE_Int iblk_lower = my_pid * loc_nblks;
    if (glb_nblks > loc_nblks * num_procs) {
        HYPRE_Int remain_nblks = glb_nblks - loc_nblks * num_procs;
        if (my_pid < remain_nblks) {
            loc_nblks++;
        }
        iblk_lower += MIN(my_pid, remain_nblks);
    }
    loc_nrows = loc_nblks * ndof;
    ilower = iblk_lower * ndof;
    iupper = ilower + loc_nrows - 1;

    if (my_pid == 0) printf("HYPRE_Int %d bytes, HYPRE_BigInt %d bytes HYPRE_Real %d bytes\n", 
        sizeof(HYPRE_Int), sizeof(HYPRE_BigInt), sizeof(HYPRE_Real));
    int ret;
    const char * dir_name = "/storage/hpcauser/zongyi/HUAWEI/baseline-hypre/data";
    char pathname[100], filename[100];
    sprintf(pathname, "%s/%s/%dx%dx%d", dir_name, case_name, gx, gy, gz);
    if (my_pid == 0) printf("reading data from %s\n", pathname);

    // 读入二进制的矩阵数据
    HYPRE_Int * dist_row_ptr = (HYPRE_Int *) malloc(sizeof(HYPRE_Int) * (loc_nrows+1));// 分布式存储的行指针（数值为全局号）

    sprintf(filename, "%s/Ai.bin", pathname);
    if ((ret = read_binary(dist_row_ptr, filename, sizeof(HYPRE_Int), ilower, loc_nrows+1)) != loc_nrows+1) {
        printf("Error! not enough rows\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    const HYPRE_Int loc_nnz = dist_row_ptr[loc_nrows] - dist_row_ptr[0];// 本进程负责的行内一共的非零元数


    // printf(" proc %d il %d iu %d nnz %d start %d\n", my_pid, ilower, iupper, loc_nnz, dist_row_ptr[0]);

    HYPRE_Int * dist_col_idx = (HYPRE_Int *) malloc(loc_nnz * sizeof(HYPRE_Int));// 分布式存储的列序号（数值为全局号）
    HYPRE_Real* dist_vals = (HYPRE_Real*)malloc(loc_nnz * sizeof(HYPRE_Real));
    sprintf(filename, "%s/Aj.bin", pathname);
    if ((ret = read_binary(dist_col_idx, filename, sizeof(HYPRE_Int), dist_row_ptr[0], loc_nnz)) != loc_nnz) {
        printf("Error! not enough dist_col_idx: %d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    sprintf(filename, "%s/Av.bin", pathname);
    if ((ret = read_binary(dist_vals, filename, sizeof(HYPRE_Real), dist_row_ptr[0], loc_nnz)) != loc_nnz) {
        printf("Error! not enough dist_vals: %d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, 3);
    }
    
    check_input(dist_vals   , loc_nnz);

    // 读入二进制的向量数据
    HYPRE_Real* dist_b = (HYPRE_Real*) malloc (loc_nrows * sizeof(HYPRE_Real));
    HYPRE_Real* dist_x = (HYPRE_Real*) malloc (loc_nrows * sizeof(HYPRE_Real));
    sprintf(filename, "%s/b.bin", pathname);
    if ((ret = read_binary(dist_b, filename, sizeof(HYPRE_Real), ilower, loc_nrows)) != loc_nrows) {
        printf("Error! not enough b\n");
        MPI_Abort(MPI_COMM_WORLD, 4);
    }

    // printf("dist_b : %lf %lf ... %lf %lf\n", dist_b[0], dist_b[1], dist_b[loc_nrows-2], dist_b[loc_nrows-1]);


    sprintf(filename, "%s/x0.bin", pathname);
    if (strstr(case_name, "FIM"  ) || (ret = read_binary(dist_x, filename, sizeof(HYPRE_Real), ilower, loc_nrows)) != loc_nrows) {
        if (my_pid == 0) printf(" ====> use zero inital guess instead!\n");
        for (HYPRE_Int i = 0; i < loc_nrows; i++)
            dist_x[i] = 0.0;
    }

    check_input(dist_b, loc_nrows);
    check_input(dist_x, loc_nrows);
    {// 自己做点积
        double self_dot[2] = {0.0, 0.0};
        for (HYPRE_BigInt i = 0; i < loc_nrows; i++) {
            self_dot[0] += (double) dist_b[i] * (double) dist_b[i];
            self_dot[1] += (double) dist_x[i] * (double) dist_x[i];
        }
        double glb_dot[2];
        MPI_Allreduce(self_dot, glb_dot, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (my_pid == 0) {
            printf("My calc dot\n");
            printf("(  b,   b) = %.15e\n",  glb_dot[0]);
            printf("(  x,   x) = %.15e\n",  glb_dot[1]);
        }
    }

    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_IJVector b, x, y;
    HYPRE_ParVector par_b, par_x, par_y;
    HYPRE_Solver solver, precond;

    // HYPRE_Init(argc, argv);
    HYPRE_Init();
    for (HYPRE_Int test = 0; test < TEST_CNT; test++) {
    
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);// Create the matrix.
    // Note that this is a square matrix, so we indicate the row partition size twice (since number of rows = number of cols)
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);// Choose a parallel csr format storage (see the User's Manual)
    HYPRE_IJMatrixInitialize(A);// Initialize before setting coefficients

    /* Create the rhs and solution */
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b); HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR); HYPRE_IJVectorInitialize(b); 
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x); HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR); HYPRE_IJVectorInitialize(x); 
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &y); HYPRE_IJVectorSetObjectType(y, HYPRE_PARCSR); HYPRE_IJVectorInitialize(y);

    {// 将矩阵和向量的数据塞给hypre
        HYPRE_Int * row_idx = (HYPRE_Int *) malloc(sizeof(HYPRE_Int) * loc_nrows);// 存本进程负责的行的全局行号
        HYPRE_Int * nnz_per_row = (HYPRE_Int *) malloc(sizeof(HYPRE_Int) * loc_nrows);// 存本进程负责的行每行有多少非零元
        for (HYPRE_Int gi = ilower; gi <= iupper; gi++) {// 全局行号
            HYPRE_Int i = gi - ilower;
            row_idx[i] = gi;
            nnz_per_row[i] = dist_row_ptr[i+1] - dist_row_ptr[i];
        }
        HYPRE_IJMatrixSetOMPFlag(A, 1);
        HYPRE_IJMatrixSetValues(A, loc_nrows, nnz_per_row, row_idx, dist_col_idx, dist_vals);

        HYPRE_IJVectorSetValues(b, loc_nrows, row_idx, dist_b);
        HYPRE_IJVectorSetValues(x, loc_nrows, row_idx, dist_x);
        free(row_idx);
        free(nnz_per_row);
    }
    
    HYPRE_IJMatrixAssemble(A);// Assemble after setting the coefficients
    HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);// Get the parcsr matrix object to use

    HYPRE_IJVectorAssemble(b);  HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorAssemble(x);  HYPRE_IJVectorGetObject(x, (void **) &par_x); 
    HYPRE_IJVectorAssemble(y);  HYPRE_IJVectorGetObject(y, (void **) &par_y);
    
    HYPRE_Real b_dot;
    HYPRE_ParVectorInnerProd(par_b, par_b, & b_dot);
// #ifdef DEBUG
    HYPRE_Real x_dot, Ab_dot, Ax_dot;
    HYPRE_ParVectorInnerProd(par_x, par_x, & x_dot);
    if (my_pid == 0) printf("(  b,   b) = %.20e\n",  b_dot);
    if (my_pid == 0) printf("(  x,   x) = %.20e\n",  x_dot);
    
    HYPRE_ParCSRMatrixMatvec(1.0, parcsr_A, par_b, 0.0, par_y);
    HYPRE_ParVectorInnerProd(par_y, par_y, &Ab_dot);
    if (my_pid == 0) printf("(A*b, A*b) = %.20e\n", Ab_dot);
    HYPRE_ParCSRMatrixMatvec(1.0, parcsr_A, par_x, 0.0, par_y);
    HYPRE_ParVectorInnerProd(par_y, par_y, &Ax_dot);
    if (my_pid == 0) printf("(A*x, A*x) = %.20e\n", Ax_dot);
// #endif

    HYPRE_Int num_iterations;
    HYPRE_Real final_res_norm;
    if (strstr(case_name, "LASER")) {/* PCG */
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
        HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
        HYPRE_PCGSetTol(solver, 1e-9); /* conv. tolerance */
        HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
        HYPRE_PCGSetPrintLevel(solver, 2); /* prints out the iteration info */
        HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

        {// 预条件AMG
            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetPrintLevel(precond, 1);
            // HYPRE_BoomerAMGSetInterpType(precond, 0);
            HYPRE_BoomerAMGSetInterpType(precond, 6);// 似乎选择默认的type 6会更好???
            HYPRE_BoomerAMGSetAggNumLevels(precond, 1);
            HYPRE_BoomerAMGSetCoarsenType(precond, 8);
            HYPRE_BoomerAMGSetMaxRowSum(precond, 0.9);
            HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25);
            HYPRE_BoomerAMGSetMaxLevels(precond, 25);
            HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
            HYPRE_BoomerAMGSetNumSweeps(precond, 1);
            HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
            HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
            HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                        (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
        }

        /* Now setup and solve! */
        HYPRE_Real t_setup = MPI_Wtime();
        HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x); t_setup = MPI_Wtime() - t_setup;
        HYPRE_Real t_calc = MPI_Wtime();
        HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x); t_calc = MPI_Wtime() - t_calc;

        /* Run info - needed logging turned on */
        HYPRE_PCGGetNumIterations(solver, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (my_pid == 0) {
            printf("\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Time cost %.5f %.5f %.5f %d\n", t_setup, t_calc, t_setup + t_calc, num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }
        /* Destroy solver and preconditioner */
        HYPRE_ParCSRPCGDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
    }
    else if (strstr(case_name, "SOLID")) {
        /* Create solver */
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
        HYPRE_PCGSetTol(solver, 1e-9); /* conv. tolerance */
        HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
        HYPRE_PCGSetPrintLevel(solver, 2); /* prints out the iteration info */
        HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

        {// 预条件AMG
            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetPrintLevel(precond, 1);
            HYPRE_BoomerAMGSetMaxLevels(precond, 25);
            HYPRE_BoomerAMGSetNumFunctions(precond, ndof);
            HYPRE_BoomerAMGSetNodal(precond, 4);// 0、1、2都不能收敛，会因为预条件的不正定而CG直接断掉，3大约125次，4、5、6均大约70次
            // HYPRE_BoomerAMGSetCoarsenType(precond, 6);// PMIS和HMIS差别不大，其它的都不收敛直接断掉
            HYPRE_BoomerAMGSetRelaxType(precond, 26);// vector version of symGS
            // HYPRE_BoomerAMGSetAggNumLevels(precond, 1);// 会报错
            // HYPRE_BoomerAMGSetRelaxOrder(precond, 1);// 用不用它对迭代数差别不大，但不用它单步时间快一些
            HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
            HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
            HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                        (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
        }
        HYPRE_Real t_setup = MPI_Wtime();
        HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x); t_setup = MPI_Wtime() - t_setup;
        HYPRE_Real t_calc = MPI_Wtime();
        HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x); t_calc = MPI_Wtime() - t_calc;

        HYPRE_PCGGetNumIterations(solver, &num_iterations);// Run info - needed logging turned on
        HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (my_pid == 0) {
            printf("\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Time cost %.5f %.5f %.5f %d\n", t_setup, t_calc, t_setup + t_calc, num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }
        HYPRE_ParCSRPCGDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
    }
    else if (strstr(case_name, "OIL")) {
        HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
        HYPRE_GMRESSetMaxIter(solver, 1000);
        HYPRE_GMRESSetPrintLevel(solver, 2); /* prints out the iteration info */
        HYPRE_GMRESSetLogging(solver, 1); /* needed to get run info later */

        if (strstr(case_name, "IMPEC")) {
            HYPRE_GMRESSetKDim(solver, 10);
            HYPRE_GMRESSetTol(solver, 1e-7);

            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetPrintLevel(precond, 1);
            HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */

            // Best Practice AMG combined with csr.fasp
            HYPRE_BoomerAMGSetNumSweeps(precond, 1);
            HYPRE_BoomerAMGSetCoarsenType(precond, 10);
            HYPRE_BoomerAMGSetAggNumLevels(precond, 1);
            HYPRE_BoomerAMGSetInterpType(precond, 6);
            HYPRE_BoomerAMGSetPMaxElmts(precond, 5);
            HYPRE_BoomerAMGSetRelaxWt(precond, 1.1);// 光滑权重
            HYPRE_BoomerAMGSetStrongThreshold(precond, 0.3);// threshold according to csr.fasp
            HYPRE_BoomerAMGSetRelaxOrder(precond, 1);// CF order according to csr.fasp
            HYPRE_BoomerAMGSetMaxRowSum(precond, 0.9);

            HYPRE_BoomerAMGSetMaxLevels(precond, 25);
            HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
            HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
            /* Set the preconditioner */
            HYPRE_GMRESSetPrecond(solver,   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
        }
        else if (strstr(case_name, "FIM"  )) {
            HYPRE_GMRESSetTol(solver, 1e-4);
            HYPRE_GMRESSetKDim(solver, 30);

            // system形式的AMG可以！！！效果也可以
            // 最好参数：nodal范数类型2或4，粗化类型6，松弛类型26（symGS，要改掉源文件不要用GE！！否则易报错）次序CF
            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetPrintLevel(precond, 1);
            HYPRE_BoomerAMGSetMaxLevels(precond, 25);
            HYPRE_BoomerAMGSetNumFunctions(precond, ndof);
            HYPRE_BoomerAMGSetNodal(precond, 2);// 0和1和2和4和5: 30，3: 52次
            HYPRE_BoomerAMGSetCoarsenType(precond, 6);// 8(PMIS) 10(HMIS)都差别不大~30次，3要40次，0和6最好都要~16次
            HYPRE_BoomerAMGSetRelaxType(precond, 26);
            HYPRE_BoomerAMGSetRelaxOrder(precond, 1);// 这个很重要
            // HYPRE_BoomerAMGSetInterpType(precond, 3);
            // HYPRE_BoomerAMGSetAggNumLevels(precond, 1);// 用这个也会报错！
            // HYPRE_BoomerAMGSetPMaxElmts(precond, 5);
            // HYPRE_BoomerAMGSetRelaxWt(precond, 1.1);// 光滑权重

            HYPRE_BoomerAMGSetTol(precond, 0.0);
            HYPRE_BoomerAMGSetMaxIter(precond, 1);
            HYPRE_GMRESSetPrecond(solver,   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
        } else assert(0);
        
        /* Now setup and solve! */
        HYPRE_Real t_setup = MPI_Wtime();
        HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x); t_setup = MPI_Wtime() - t_setup;
        HYPRE_Real t_calc = MPI_Wtime();
        HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_b, par_x); t_calc = MPI_Wtime() - t_calc;

        /* Run info - needed logging turned on */
        HYPRE_GMRESGetNumIterations(solver, &num_iterations);
        HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (my_pid == 0) {
            printf("\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Time cost %.5f %.5f %.5f %d\n", t_setup, t_calc, t_setup + t_calc, num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }
        
        HYPRE_ParCSRGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
    }
    else if (strstr(case_name, "GRAPES")) {
        HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
        HYPRE_GMRESSetKDim(solver, 5);
        HYPRE_GMRESSetMaxIter(solver, 1000);
        HYPRE_GMRESSetPrintLevel(solver, 2); /* prints out the iteration info */
        HYPRE_GMRESSetLogging(solver, 1); /* needed to get run info later */
        HYPRE_GMRESSetTol(solver, 0.0);
        HYPRE_GMRESSetAbsoluteTol(solver, 1e-5);
        {// 预条件AMG
            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetPrintLevel(precond, 1);
            HYPRE_BoomerAMGSetMaxLevels(precond, 25);
            
            HYPRE_BoomerAMGSetInterpType(precond, 6);// BEST: extended classical modiﬁed interpolation (interpolation type 6)
            HYPRE_BoomerAMGSetCoarsenType(precond, 10);// BEST: HMIS coarsening (coarsen type 10),
            HYPRE_BoomerAMGSetAggNumLevels(precond, 1);// BEST: aggressive coarsening on the ﬁnest level (aggressive num levels 1)
            HYPRE_BoomerAMGSetPMaxElmts(precond, 5);// BEST: truncated interpolation to ﬁve entries per row (PMaxElmts option equals 5)
            HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25);// BEST: strength-of-connection tolerance of 0.25
            HYPRE_BoomerAMGSetNumSweeps(precond, 1);
            HYPRE_BoomerAMGSetRelaxType(precond, 6);// BEST: hybrid symmetric Gauss–Seidel (relax type 6)

            HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
            HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
            HYPRE_GMRESSetPrecond(solver,   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
        }
        HYPRE_Real t_setup = MPI_Wtime();
        HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x); t_setup = MPI_Wtime() - t_setup;
        HYPRE_Real t_calc = MPI_Wtime();
        HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_b, par_x); t_calc = MPI_Wtime() - t_calc;

        /* Run info - needed logging turned on */
        HYPRE_GMRESGetNumIterations(solver, &num_iterations);
        HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (my_pid == 0) {
            printf("\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Time cost %.5f %.5f %.5f %d\n", t_setup, t_calc, t_setup + t_calc, num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }
        
        HYPRE_ParCSRGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
    }
    else assert(0);

    // 计算真实残差
    HYPRE_ParVectorCopy(par_b, par_y);// y = b
    HYPRE_ParCSRMatrixMatvec(-1.0, parcsr_A, par_x, 1.0, par_y);// y += -A*x
    HYPRE_Real r_dot;
    HYPRE_ParVectorInnerProd(par_y, par_y, &r_dot);
    if (my_pid == 0) printf("\033[1;35mtrue ||r||/||b||= %20.16e\033[0m\n", sqrt(r_dot/b_dot));

    /* Clean up */
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b); HYPRE_IJVectorDestroy(x); HYPRE_IJVectorDestroy(y);
    
    }// test loop

    HYPRE_Finalize();
    MPI_Finalize();

    free(dist_b); free(dist_x);
    free(dist_row_ptr); free(dist_col_idx); free(dist_vals);

    return 0;
}