#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"
// #include "HYPRE_FEI.h"
// #include "HYPRE_LSI_mli.h"
// #include "HYPRE_"

#define data_t double
#define MIN(a, b) (a) < (b) ? (a) : (b)
#ifndef TEST_CNT
#define TEST_CNT 5
#endif

int read_binary(void * buf, const char * filename, int size_of_elems, long start, int nums) {
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
    int ret;
    
    ret = fread(buf, size_of_elems, nums, fp);
    return ret;
}

void check_input(const data_t * input, int num) {
    for (int i = 0; i < num; i++)
        // assert(input[i] == input[i]);
        assert(abs(input[i]) < 1e20);
}

int main(int argc, char * argv[])
{
    int my_pid, num_procs;
    int glb_nrows;
    int glb_nblks = -1, ndof = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);

    int cnt = 1;
    const char * case_name = argv[cnt++];
    int gx = atoi(argv[cnt++]),
        gy = atoi(argv[cnt++]),
        gz = atoi(argv[cnt++]);
    if (strcmp(case_name, "LASER") == 0) {
        glb_nrows = gx * gy * gz;
    }
    else if (strcmp(case_name, "LASER-3T") == 0) {
        ndof = 3;
        glb_nblks = gx * gy * gz;
        glb_nrows = glb_nblks * ndof;
    }
    else if (strcmp(case_name, "OIL-IMPEC") == 0) {
        glb_nrows =  gx * gy * gz + 5;
    }
    else if (strcmp(case_name, "OIL-FIM") == 0) {
        int num_wells = 5;
        if (my_pid == 0) printf("OIL-FIM num wells: %d\n", num_wells);
        glb_nblks = gx * gy * gz + num_wells;
        ndof = 4;
        glb_nrows = glb_nblks * ndof;
    }
    else {
        if (my_pid == 0) printf("Error: Unrecognized case name %s\n", case_name);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }


    int loc_nrows;
    int loc_nblks = -1;
    int ilower, iupper;// 本进程负责的起始行号 和末尾行号（闭区间）
    if (glb_nblks != -1) {// 此时是多自由度: OIL-FIM，需要先对blk做划分
        if (my_pid == 0) printf("using blk-partition!\n");
        loc_nblks = glb_nblks / num_procs;
        int iblk_lower = my_pid * loc_nblks;
        if (glb_nblks > loc_nblks * num_procs) {
            int remain_nblks = glb_nblks - loc_nblks * num_procs;
            if (my_pid < remain_nblks) {
                loc_nblks++;
            }
            iblk_lower += MIN(my_pid, remain_nblks);
        }
        loc_nrows = loc_nblks * ndof;
        ilower = iblk_lower * ndof;
        iupper = ilower + loc_nrows - 1;
        // printf("proc %d blk [%d,%d), ilower %d iupper %d\n", my_pid, iblk_lower, iblk_lower+loc_nblks, ilower, iupper);
    } 
    else {// 直接对矩阵的各行做划分
        loc_nrows = glb_nrows / num_procs;// 本进程负责的行数
        ilower = my_pid * loc_nrows;
        if (glb_nrows > loc_nrows * num_procs) {// 各自再分一些行
            int remain_nrows = glb_nrows - loc_nrows * num_procs;// 余下没分的行数
            if (my_pid < remain_nrows) {
                loc_nrows ++;
            }
            ilower += MIN(my_pid, remain_nrows);
        }
        iupper = ilower + loc_nrows - 1;// 本进程负责的末尾行号（inclusive，闭区间）
    }

    int ret;
    const char * dir_name = "/storage/hpcauser/zongyi/HUAWEI/baseline-hypre/data";
    char pathname[100], filename[100];
    sprintf(pathname, "%s/%s/%dx%dx%d", dir_name, case_name, gx, gy, gz);
    if (my_pid == 0) printf("reading data from %s\n", pathname);

    // 读入二进制的矩阵数据
    int * dist_row_ptr = (int *) malloc(sizeof(int) * (loc_nrows+1));// 分布式存储的行指针（数值为全局号）

    sprintf(filename, "%s/Ai.bin", pathname);
    if ((ret = read_binary(dist_row_ptr, filename, sizeof(int), ilower, loc_nrows+1)) != loc_nrows+1) {
        printf("Error! not enough rows\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    const int loc_nnz = dist_row_ptr[loc_nrows] - dist_row_ptr[0];// 本进程负责的行内一共的非零元数


    // printf(" proc %d il %d iu %d nnz %d start %d\n", my_pid, ilower, iupper, loc_nnz, dist_row_ptr[0]);

    int * dist_col_idx = (int *) malloc(loc_nnz * sizeof(int));// 分布式存储的列序号（数值为全局号）
    data_t* dist_vals = (data_t*)malloc(loc_nnz * sizeof(data_t));
    sprintf(filename, "%s/Aj.bin", pathname);
    if ((ret = read_binary(dist_col_idx, filename, sizeof(int), dist_row_ptr[0], loc_nnz)) != loc_nnz) {
        printf("Error! not enough dist_col_idx: %d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    // assert(num_procs == 1);
    // FILE * fp = fopen("colidx", "w+");
    // for (int i = 0; i < loc_nnz; i++)
    //     fprintf(fp, "%d\n", dist_col_idx[i]);
    // fclose(fp);

    sprintf(filename, "%s/Av.bin", pathname);
    if ((ret = read_binary(dist_vals, filename, sizeof(data_t), dist_row_ptr[0], loc_nnz)) != loc_nnz) {
        printf("Error! not enough dist_vals: %d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, 3);
    }
    
    check_input(dist_vals   , loc_nnz);

    // 读入二进制的向量数据
    data_t* dist_b = (data_t*) malloc (loc_nrows * sizeof(data_t));
    data_t* dist_x = (data_t*) malloc (loc_nrows * sizeof(data_t));
    sprintf(filename, "%s/b.bin", pathname);
    if ((ret = read_binary(dist_b, filename, sizeof(data_t), ilower, loc_nrows)) != loc_nrows) {
        printf("Error! not enough b\n");
        MPI_Abort(MPI_COMM_WORLD, 4);
    }

    // printf("dist_b : %lf %lf ... %lf %lf\n", dist_b[0], dist_b[1], dist_b[loc_nrows-2], dist_b[loc_nrows-1]);


    sprintf(filename, "%s/x0.bin", pathname);
    if (strstr(case_name, "FIM"  ) || (ret = read_binary(dist_x, filename, sizeof(data_t), ilower, loc_nrows)) != loc_nrows) {
        if (my_pid == 0) printf(" ====> use zero inital guess instead!\n");
        for (int i = 0; i < loc_nrows; i++)
            dist_x[i] = 0.0;
    }

    check_input(dist_b, loc_nrows);
    check_input(dist_x, loc_nrows);

    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_IJVector b, x, y;
    HYPRE_ParVector par_b, par_x, par_y;
    HYPRE_Solver solver, precond;

    // HYPRE_Init(argc, argv);
    HYPRE_Init();
    for (int test = 0; test < TEST_CNT; test++) {
    
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);// Create the matrix.
    // Note that this is a square matrix, so we indicate the row partition size twice (since number of rows = number of cols)
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);// Choose a parallel csr format storage (see the User's Manual)
    HYPRE_IJMatrixInitialize(A);// Initialize before setting coefficients

    /* Create the rhs and solution */
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b); HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR); HYPRE_IJVectorInitialize(b); 
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x); HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR); HYPRE_IJVectorInitialize(x); 
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &y); HYPRE_IJVectorSetObjectType(y, HYPRE_PARCSR); HYPRE_IJVectorInitialize(y);

    {// 将矩阵和向量的数据塞给hypre
        int * row_idx = (int *) malloc(sizeof(int) * loc_nrows);// 存本进程负责的行的全局行号
        int * nnz_per_row = (int *) malloc(sizeof(int) * loc_nrows);// 存本进程负责的行每行有多少非零元
        for (int gi = ilower; gi <= iupper; gi++) {// 全局行号
            int i = gi - ilower;
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
    
    double b_dot;
    HYPRE_ParVectorInnerProd(par_b, par_b, & b_dot);
// #ifdef DEBUG
    double x_dot, Ab_dot;
    HYPRE_ParVectorInnerProd(par_x, par_x, & x_dot);
    if (my_pid == 0) printf("(  b,   b) = %.20e\n",  b_dot);
    if (my_pid == 0) printf("(  x,   x) = %.20e\n",  x_dot);
    
    HYPRE_ParCSRMatrixMatvec(1.0, parcsr_A, par_b, 0.0, par_y);
    HYPRE_ParVectorInnerProd(par_y, par_y, &Ab_dot);
    if (my_pid == 0) printf("(A*b, A*b) = %.20e\n", Ab_dot);

// #endif

    if (strstr(case_name, "LASER")) {/* PCG */
        int num_iterations;
        double final_res_norm;
        /* Create solver */
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
        HYPRE_PCGSetTol(solver, 1e-9); /* conv. tolerance */
        HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
        HYPRE_PCGSetPrintLevel(solver, 2); /* prints out the iteration info */
        HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

        // if (strstr(case_name, "3T")) {
        //     HYPRE_BoomerAMGCreate(&precond);
        //     HYPRE_BoomerAMGSetPrintLevel(precond, 1);
        //     HYPRE_BoomerAMGSetNumFunctions(precond, ndof);
        //     // HYPRE_BoomerAMGSetNodal(precond, 2);
        //     // HYPRE_BoomerAMGSetCoarsenType(precond, 6);// 8(PMIS) 10(HMIS)都差别不大~30次，3要40次，0和6最好都要~16次
        //     HYPRE_BoomerAMGSetInterpType(precond, 6);
        //     HYPRE_BoomerAMGSetRelaxType(precond, 26);
        //     HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        //     HYPRE_BoomerAMGSetTol(precond, 0.0);
        //     HYPRE_BoomerAMGSetMaxIter(precond, 1);
        //     HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
        //                                 (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
        // } else
        {// 预条件AMG
            /* Now set up the AMG preconditioner and specify any parameters */
            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
            // HYPRE_BoomerAMGSetInterpType(precond, 0);
            HYPRE_BoomerAMGSetInterpType(precond, 6);// 似乎选择默认的type 6会更好???
            // HYPRE_BoomerAMGSetAggNumLevels(precond, 1);
            HYPRE_BoomerAMGSetCoarsenType(precond, 8);
            HYPRE_BoomerAMGSetMaxRowSum(precond, 0.9);
            HYPRE_BoomerAMGSetStrongThreshold(precond, 0.25);
            HYPRE_BoomerAMGSetMaxLevels(precond, 10);
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
    else if (strstr(case_name, "OIL")) {
        int num_iterations;
        double final_res_norm;

        HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
        HYPRE_GMRESSetMaxIter(solver, 1000);
        HYPRE_GMRESSetPrintLevel(solver, 2); /* prints out the iteration info */
        HYPRE_GMRESSetLogging(solver, 1); /* needed to get run info later */

        if (strstr(case_name, "IMPEC")) {
            HYPRE_GMRESSetKDim(solver, 10);
            HYPRE_GMRESSetTol(solver, 1e-7);

            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
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

            HYPRE_BoomerAMGSetMaxLevels(precond, 20);
            HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
            HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
            /* Set the preconditioner */
            HYPRE_GMRESSetPrecond(solver,   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
        }
        else if (strstr(case_name, "FIM"  )) {
            HYPRE_GMRESSetTol(solver, 1e-4);
            HYPRE_GMRESSetKDim(solver, 30);

            /* // 搞不定
            HYPRE_LSI_MLICreate(MPI_COMM_WORLD, &precond);
            HYPRE_LSI_MLISetParams(precond, "MLI outputLevel 2");
            HYPRE_LSI_MLIGetNullSpace
            // HYPRE_LSI_MLISetParams(precond, "MLI numLevels 20");
            // HYPRE_LSI_MLISetParams(precond, "MLI maxIterations 1");
            // HYPRE_LSI_MLISetParams(precond, "MLI cycleType V");
            HYPRE_LSI_MLISetParams(precond, "MLI strengthThreshold 0.08");
            // HYPRE_LSI_MLISetParams(precond, "MLI method AMGSAe");
            // HYPRE_LSI_MLISetParams(precond, "MLI coarsenScheme falgout");
            // HYPRE_LSI_MLISetParams(precond, "MLI smoother GS");
            // HYPRE_LSI_MLISetParams(precond, "MLI coarseSolver SuperLU");
            // HYPRE_LSI_MLISetParams(precond, "MLI numSweeps 2");
            HYPRE_LSI_MLISetParams(precond, "MLI nullSpaceDim 100");
            HYPRE_LSI_MLISetParams(precond, "MLI printNullSpace");
            // HYPRE_LSI_MLISetParams(precond, "MLI Pweight 1.2");
            // HYPRE_LSI_MLISetParams(precond, "MLI smootherWeight 5");
            HYPRE_GMRESSetPrecond(solver,   (HYPRE_PtrToSolverFcn) HYPRE_LSI_MLISolve,
                                            (HYPRE_PtrToSolverFcn) HYPRE_LSI_MLISetup, precond);
            */

            /*
            HYPRE_MGRCreate(&precond);
            HYPRE_MGRSetPrintLevel(precond, 1);
            HYPRE_MGRSetBlockSize(precond, 4);
            HYPRE_MGRSetCoarseGridMethod(precond, 0);
            HYPRE_MGRSetCoarseGridPrintLevel(precond, 2);
            HYPRE_MGRSetTol(precond, 0.0);
            HYPRE_MGRSetMaxIter(precond, 1);
            HYPRE_MGRSetCpointsByPointMarkerArray
            */

            // system形式的AMG可以！！！效果也可以
            // 最好参数：nodal范数类型2或4，粗化类型6，松弛类型26（symGS，要改掉源文件不要用GE！！否则易报错）次序CF
            HYPRE_BoomerAMGCreate(&precond);
            HYPRE_BoomerAMGSetPrintLevel(precond, 1);
            HYPRE_BoomerAMGSetNumFunctions(precond, ndof);
            HYPRE_BoomerAMGSetNodal(precond, 2);// 0和1和2和4和5: 30，3: 52次
            HYPRE_BoomerAMGSetCoarsenType(precond, 6);// 8(PMIS) 10(HMIS)都差别不大~30次，3要40次，0和6最好都要~16次
            HYPRE_BoomerAMGSetRelaxType(precond, 26);
            HYPRE_BoomerAMGSetRelaxOrder(precond, 1);// 这个很重要
            // HYPRE_BoomerAMGSetInterpType(precond, 3);
            // HYPRE_BoomerAMGSetAggNumLevels(precond, 1);// 用这个也会报错！
            // HYPRE_BoomerAMGSetPMaxElmts(precond, 5);
            // HYPRE_BoomerAMGSetRelaxWt(precond, 1.1);// 光滑权重

            // ILU的不行
            // HYPRE_BoomerAMGSetSmoothNumLevels(precond, 6);
            // HYPRE_BoomerAMGSetSmoothType(precond, 5);
            // HYPRE_BoomerAMGSetSmoothNumSweeps(precond, 1);
            // HYPRE_ILUSetType(precond, 30);
            // HYPRE_ILUSetLevelOfFill(precond, 1);
            // HYPRE_BoomerAMGSetILUType(precond, 30);
            // HYPRE_BoomerAMGSetILULevel(precond, 1);
            // HYPRE_BoomerAMGSetILUMaxIter(precond, 1);

            HYPRE_BoomerAMGSetTol(precond, 0.0);
            HYPRE_BoomerAMGSetMaxIter(precond, 1);
            HYPRE_GMRESSetPrecond(solver,   (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                                            (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
            
            /*
            //可以，配合重启长度为30更佳，但也只能对第0时刻的算例有效
            HYPRE_ILUCreate(&precond);
            HYPRE_ILUSetPrintLevel(precond, 1);
            HYPRE_ILUSetTol(precond, 0.0);
            HYPRE_ILUSetMaxIter(precond, 1);
            HYPRE_ILUSetType(precond, 30);
            HYPRE_ILUSetLevelOfFill(precond, 2);
            HYPRE_GMRESSetPrecond(solver,   (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve,
                                            (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup, precond);
            */
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

        if (strstr(case_name, "IMPEC")) HYPRE_BoomerAMGDestroy(precond);
        if (strstr(case_name, "FIM"  )) HYPRE_BoomerAMGDestroy(precond);
        // if (strstr(case_name, "FIM"))   HYPRE_LSI_MLIDestroy(precond);
    }
    else assert(0);

    // 计算真实残差
    HYPRE_ParVectorCopy(par_b, par_y);// y = b
    HYPRE_ParCSRMatrixMatvec(-1.0, parcsr_A, par_x, 1.0, par_y);// y += -A*x
    double r_dot;
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