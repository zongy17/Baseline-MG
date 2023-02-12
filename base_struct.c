#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "HYPRE_struct_ls.h"
#include "_hypre_struct_ls.h"

#define MIN(a, b) (a) < (b) ? (a) : (b)
#ifndef TEST_CNT
#define TEST_CNT 1
#endif

int main(int argc, char * argv[])
{
    int my_pid, num_procs;
    HYPRE_Int glb_nrows;

    MPI_Init(&argc, &argv);

    HYPRE_Int cnt = 1;
    const char * case_name = argv[cnt++];
    HYPRE_Int gx = atoi(argv[cnt++]),
        gy = atoi(argv[cnt++]),
        gz = atoi(argv[cnt++]);
    HYPRE_Int px = atoi(argv[cnt++]),
        py = atoi(argv[cnt++]),
        pz = atoi(argv[cnt++]);
    const char * precond_name = argv[cnt++];
    if (strcmp(case_name, "LASER") == 0) {
        glb_nrows = gx * gy * gz;
    }
    else if (strcmp(case_name, "GRAPES") == 0) {
        assert(pz == 1);
        glb_nrows = gx * gy * gz;
    }
    else {
        MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
        if (my_pid == 0) printf("Error: Unrecognized case name %s\n", case_name);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    HYPRE_Int lx, ly, lz;
    if (gx % px != 0 || gy % py != 0 || gz % pz != 0) {// check if input is legal
        MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
        if (my_pid == 0) printf("Error: Illegal input! gx/y/z must be fully divided by px/y/z!\n", gx, gy, gz, px, py, pz);
        MPI_Abort(MPI_COMM_WORLD, -1);
    } else {
        lx = gx / px;
        ly = gy / py;
        lz = gz / pz;
    }
    const HYPRE_Int ndim = 3;

    MPI_Comm cart_comm;
    int cart_ids[ndim], pidx, pidy, pidz;
    HYPRE_Int ilower[ndim], iupper[ndim];// 注意i[2]在最外维，i[0]在最内维 所以令0: z, 1: x, 2: y
    {// Setup Cartesian comm
        HYPRE_Int dims[ndim] = {py, px, pz};// 注意这里保持一致！！！
        HYPRE_Int periods[ndim] = {0, 0, 0};
        MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, 0, &cart_comm);
        MPI_Comm_rank(cart_comm, &my_pid);
        MPI_Comm_size(cart_comm, &num_procs);
        MPI_Cart_coords(cart_comm, my_pid, ndim, cart_ids);
        pidx = cart_ids[1];
        pidy = cart_ids[0];
        pidz = cart_ids[2];

        ilower[0] = pidz * lz; iupper[0] = ilower[0] + lz - 1;
        ilower[1] = pidx * lx; iupper[1] = ilower[1] + lx - 1;
        ilower[2] = pidy * ly; iupper[2] = ilower[2] + ly - 1;
    }
    HYPRE_Init();
    HYPRE_StructGrid     grid;
    HYPRE_StructStencil  stencil;
    HYPRE_StructMatrix   A;
    HYPRE_StructVector   b, x, y;
    HYPRE_StructSolver   solver;
    HYPRE_StructSolver   precond;

    for (HYPRE_Int test = 0; test < TEST_CNT; test++) {
    
    HYPRE_Int num_diag = -1;
    if (strcmp(case_name, "LASER") == 0) {
        assert(sizeof(HYPRE_Real) == sizeof(double));
        {// Setup a grid
            HYPRE_StructGridCreate(cart_comm, ndim, &grid);// Create an empty 3D grid object
            HYPRE_StructGridSetExtents(grid, ilower, iupper);// Add a new box to the grid
            HYPRE_StructGridAssemble(grid);// a collective call finalizing the grid assembly.
        }
        num_diag = 7;
        HYPRE_Int offsets[num_diag][ndim] = {{0,0,-1}, {0,-1,0}, {-1,0,0}, {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
        {// Define the discretization stencil
            HYPRE_StructStencilCreate(ndim, num_diag, &stencil);
            // Define the geometry of the stencil
            for (HYPRE_Int entry = 0; entry < num_diag; entry++)
                HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
        }
    }
    else if (strcmp(case_name, "GRAPES") == 0) {
        assert(sizeof(HYPRE_Real) == sizeof(float));
        {// Setup a grid
            HYPRE_StructGridCreate(cart_comm, ndim, &grid);// Create an empty 3D grid object
            HYPRE_StructGridSetExtents(grid, ilower, iupper);// Add a new box to the grid
            HYPRE_StructGridAssemble(grid);// a collective call finalizing the grid assembly.
        }
        num_diag = 19;
        HYPRE_Int offsets[num_diag][ndim] = {
                        { 0,-1,-1}, 
            {-1, 0,-1},{ 0, 0,-1}, { 1, 0,-1},
                        { 0, 1,-1},
            
            {-1,-1, 0},{ 0,-1, 0},{ 1, -1, 0},
            {-1, 0, 0},{ 0, 0, 0},{ 1,  0, 0},
            {-1, 1, 0},{ 0, 1, 0},{ 1,  1, 0},

                    { 0,-1, 1},
            {-1, 0, 1},{ 0, 0, 1},{ 1, 0, 1},
                    { 0, 1, 1}
        };
        {// Define the discretization stencil
            HYPRE_StructStencilCreate(ndim, num_diag, &stencil);
            // Define the geometry of the stencil
            for (HYPRE_Int entry = 0; entry < num_diag; entry++)
                HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
        }
    }
    
    assert(num_diag != -1);
    HYPRE_Int stencil_indices[num_diag];
    for (HYPRE_Int j = 0; j < num_diag; j++)
        stencil_indices[j] = j;
    {// Set up a Struct Matrix and Struct Vector
        HYPRE_Int tot_len = num_diag * lx * ly * lz;
        HYPRE_Real * buf = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * tot_len);
        // 首先MPI-IO读入矩阵数据
        if      (strcmp(case_name, "LASER") == 0) assert(sizeof(double) == sizeof(HYPRE_Real));
        else if (strcmp(case_name, "GRAPES")== 0) assert(sizeof(float ) == sizeof(HYPRE_Real));
        MPI_File fh = MPI_FILE_NULL;// 文件句柄
        MPI_Datatype etype = (sizeof(HYPRE_Real) == 4) ? MPI_FLOAT : MPI_DOUBLE;
        MPI_Datatype read_type = MPI_DATATYPE_NULL;// 写出类型
        HYPRE_Int size[4], subsize[4], start[4];
        size   [0] = gy       ; size   [1] = gx       ; size   [2] = gz       ; size   [3] = num_diag;
        subsize[0] = ly       ; subsize[1] = lx       ; subsize[2] = lz       ; subsize[3] = num_diag;
        start  [0] = pidy * ly; start  [1] = pidx * lx; start  [2] = pidz * lz; start  [3] = 0       ;

        MPI_Type_create_subarray(4, size, subsize, start, MPI_ORDER_C, etype, &read_type);
        MPI_Type_commit(&read_type);
        MPI_Offset displacement = 0;
        displacement *= sizeof(HYPRE_Real);// 位移要以字节为单位！
        MPI_Status status;

        int ret;
        ret = MPI_File_open(cart_comm, "mat.AOS.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if (ret != MPI_SUCCESS) {
            printf("Could not open file: ret %d\n", ret);
        }
        MPI_File_set_view(fh, displacement, etype, read_type, "native", MPI_INFO_NULL);
        MPI_File_read_all(fh, buf, tot_len, etype, &status);
        MPI_File_close(&fh);
        // 然后 传给结构化矩阵A
        HYPRE_StructMatrixCreate(cart_comm, grid, stencil, &A);// Create an empty matrix object
        HYPRE_StructMatrixInitialize(A);// Indicate that the matrix coefficients are ready to be set 
        HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, num_diag, stencil_indices, buf);
        HYPRE_StructMatrixAssemble(A);// a collective call finalizing the matrix assembly.

        MPI_Type_free(&read_type);
        free(buf);

        // 向量
        HYPRE_StructVectorCreate(cart_comm, grid, &b);// Create an empty vector object
        HYPRE_StructVectorCreate(cart_comm, grid, &x);
        HYPRE_StructVectorCreate(cart_comm, grid, &y);
        HYPRE_StructVectorInitialize(b);// Indicate that the vector coefficients are ready to be set
        HYPRE_StructVectorInitialize(x);
        HYPRE_StructVectorInitialize(y);

        tot_len = lx * ly * lz;
        buf = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * tot_len);
        MPI_Type_create_subarray(3, size, subsize, start, MPI_ORDER_C, etype, &read_type);
        MPI_Type_commit(&read_type);
        {// b
            ret = MPI_File_open(cart_comm, "b.AOS.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
            if (ret != MPI_SUCCESS) {
                printf("Could not open file: ret %d\n", ret);
            }
            MPI_File_set_view(fh, displacement, etype, read_type, "native", MPI_INFO_NULL);
            MPI_File_read_all(fh, buf, tot_len, etype, &status);
            MPI_File_close(&fh);

            HYPRE_StructVectorSetBoxValues(b, ilower, iupper, buf);
        }
        {// x
            ret = MPI_File_open(cart_comm, "x.AOS.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
            if (ret != MPI_SUCCESS) {
                printf("Could not open file: ret %d\n", ret);
            }
            MPI_File_set_view(fh, displacement, etype, read_type, "native", MPI_INFO_NULL);
            MPI_File_read_all(fh, buf, tot_len, etype, &status);
            MPI_File_close(&fh);

            HYPRE_StructVectorSetBoxValues(x, ilower, iupper, buf);
        }
        
        HYPRE_StructVectorAssemble(b);
        HYPRE_StructVectorAssemble(x);
        HYPRE_StructVectorAssemble(y);// 对y也要创建
        
        MPI_Type_free(&read_type);
        free(buf);
    }

    HYPRE_Real b_dot;
    b_dot = hypre_StructKrylovInnerProd(b, b);
    if (my_pid == 0) printf("(  b,   b) = %.20e\n",  b_dot);
    {// Check for DEBUG purpose
        HYPRE_Real x_dot, Ab_dot, Ax_dot;
        x_dot = hypre_StructKrylovInnerProd(x, x);
        if (my_pid == 0) printf("(  x,   x) = %.20e\n",  x_dot);
        HYPRE_StructMatrixMatvec(1.0, A, x, 0.0, y);
        Ax_dot= hypre_StructKrylovInnerProd(y, y);
        if (my_pid == 0) printf("(A*x, A*x) = %.20e\n", Ax_dot);
        HYPRE_StructMatrixMatvec(1.0, A, b, 0.0, y);
        Ab_dot= hypre_StructKrylovInnerProd(y, y);
        if (my_pid == 0) printf("(A*b, A*b) = %.20e\n", Ab_dot);
    }
    
    if (strcmp(case_name, "LASER") == 0) {// Setup a struct solver, preconditioner and RUN
        HYPRE_StructPCGCreate(cart_comm, &solver);
        HYPRE_StructPCGSetMaxIter(solver, 1000);
        HYPRE_StructPCGSetTol(solver, 1.0e-09);

        HYPRE_StructPCGSetTwoNorm(solver, 1 );
        HYPRE_StructPCGSetRelChange(solver, 0 );
        HYPRE_StructPCGSetPrintLevel(solver, 2);
        HYPRE_StructPCGSetLogging(solver, 1);

        if (strcmp(precond_name, "SMG") == 0) {
            if (my_pid == 0) printf("using SMG...\n");
            HYPRE_StructSMGCreate(cart_comm, &precond);
            HYPRE_StructSMGSetPrintLevel(precond, 10);
            HYPRE_StructSMGSetLogging(precond, 0);
            // HYPRE_StructSMGSetMemoryUse(precond, 0);
            HYPRE_StructSMGSetMaxIter(precond, 1);
            HYPRE_StructSMGSetTol(precond, 0.0);

            HYPRE_StructSMGSetZeroGuess(precond);
            // SMG的松弛类型是不能变的，它做的是一维半粗化，然后每个面做plane-solve，其中的plane-solve又来一遍半粗化和线松弛
            HYPRE_StructSMGSetNumPreRelax(precond, 1);
            HYPRE_StructSMGSetNumPostRelax(precond, 1);

            HYPRE_StructPCGSetPrecond(solver,   HYPRE_StructSMGSolve,
                                                HYPRE_StructSMGSetup, precond);
        }
        else if (strcmp(precond_name, "PFMG") == 0) {
            if (my_pid == 0) printf("using PFMG...\n");
            HYPRE_StructPFMGCreate(cart_comm, &precond);
            HYPRE_StructPFMGSetPrintLevel(precond, 3);
            HYPRE_StructPFMGSetLogging(precond, 0);

            HYPRE_StructPFMGSetMaxIter(precond, 1);
            HYPRE_StructPFMGSetTol(precond, 0.0);
            HYPRE_StructPFMGSetZeroGuess(precond);
            // HYPRE_StructPFMGSetRAPType(precond, 1);
            HYPRE_StructPFMGSetRelaxType(precond, 1);// 0: 655次，1: 70次（如果搭配RAPType=1则要102次），2：93次，3：89次（2和3不能使用RAPType=0）
            HYPRE_StructPFMGSetNumPreRelax(precond, 1);
            HYPRE_StructPFMGSetNumPostRelax(precond, 1);
            HYPRE_StructPCGSetPrecond(solver,   HYPRE_StructPFMGSolve,
                                                HYPRE_StructPFMGSetup, precond);
        }
        else if (strcmp(precond_name, "NOPC") == 0) { 
            if (my_pid == 0) printf("  NO precond!\n");
        }
        else {
            if (my_pid == 0) printf("  Error: Uncognized PC name: %s!\n", precond_name);
            MPI_Abort(MPI_COMM_WORLD, -101);
        }

        HYPRE_Real t_setup = MPI_Wtime();
        HYPRE_StructPCGSetup(solver, A, b, x ); t_setup = MPI_Wtime() - t_setup;
        HYPRE_Real t_calc  = MPI_Wtime();
        HYPRE_StructPCGSolve(solver, A, b, x);  t_calc  = MPI_Wtime() - t_calc;

        /* Get info and release memory */
        HYPRE_Int num_iterations;
        HYPRE_Real final_res_norm;
        HYPRE_StructPCGGetNumIterations( solver, &num_iterations );
        HYPRE_StructPCGGetFinalRelativeResidualNorm( solver, &final_res_norm );

        if (my_pid == 0) {
            printf("\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Time cost %.5f %.5f %.5f %d\n", t_setup, t_calc, t_setup + t_calc, num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }

        // Destroy solver and preconditioner
        HYPRE_StructPCGDestroy(solver);
        if      (strcmp(precond_name, "SMG" ) == 0) HYPRE_StructSMGDestroy (precond);
        else if (strcmp(precond_name, "PFMG") == 0) HYPRE_StructPFMGDestroy(precond);
    }
    else if (strcmp(case_name, "GRAPES") == 0) {
        HYPRE_StructGMRESCreate(cart_comm, &solver);
        HYPRE_StructGMRESSetKDim(solver, 10);
        HYPRE_StructGMRESSetMaxIter(solver, 1000);
        HYPRE_StructGMRESSetAbsoluteTol(solver, 1.0e-05);
        HYPRE_StructGMRESSetTol(solver, 0.0);
        HYPRE_StructGMRESSetPrintLevel(solver, 2);
        HYPRE_StructGMRESSetLogging(solver, 1);

        if (strcmp(precond_name, "SMG") == 0) {
            if (my_pid == 0) printf("using SMG...\n");
            HYPRE_StructSMGCreate(cart_comm, &precond);
            HYPRE_StructSMGSetPrintLevel(precond, 10);
            HYPRE_StructSMGSetLogging(precond, 0);
            // HYPRE_StructSMGSetMemoryUse(precond, 0);
            HYPRE_StructSMGSetMaxIter(precond, 1);
            HYPRE_StructSMGSetTol(precond, 0.0);

            HYPRE_StructSMGSetZeroGuess(precond);
            // SMG的松弛类型是不能变的，它做的是一维半粗化，然后每个面做plane-solve，其中的plane-solve又来一遍半粗化和线松弛
            HYPRE_StructSMGSetNumPreRelax(precond, 1);
            HYPRE_StructSMGSetNumPostRelax(precond, 1);

            HYPRE_StructGMRESSetPrecond(solver,   HYPRE_StructSMGSolve,
                                                HYPRE_StructSMGSetup, precond);
        }
        else if (strcmp(precond_name, "PFMG") == 0) {
            if (my_pid == 0) printf("using PFMG...\n");
            HYPRE_StructPFMGCreate(cart_comm, &precond);
            HYPRE_StructPFMGSetPrintLevel(precond, 3);
            HYPRE_StructPFMGSetLogging(precond, 0);

            HYPRE_StructPFMGSetMaxIter(precond, 1);
            HYPRE_StructPFMGSetTol(precond, 0.0);
            HYPRE_StructPFMGSetZeroGuess(precond);
            // HYPRE_StructPFMGSetRAPType(precond, 0);
            HYPRE_StructPFMGSetRelaxType(precond, 2);// 0: 655次，1: 70次，
            HYPRE_StructPFMGSetNumPreRelax(precond, 2);
            HYPRE_StructPFMGSetNumPostRelax(precond, 2);
            HYPRE_StructGMRESSetPrecond(solver,   HYPRE_StructPFMGSolve,
                                                HYPRE_StructPFMGSetup, precond);
        }
        else if (strcmp(precond_name, "NOPC") == 0) { 
            if (my_pid == 0) printf("  NO precond!\n");
        }
        else {
            if (my_pid == 0) printf("  Error: Uncognized PC name: %s!\n", precond_name);
            MPI_Abort(MPI_COMM_WORLD, -101);
        }

        HYPRE_Real t_setup = MPI_Wtime();
        HYPRE_StructGMRESSetup(solver, A, b, x ); t_setup = MPI_Wtime() - t_setup;
        HYPRE_Real t_calc  = MPI_Wtime();
        HYPRE_StructGMRESSolve(solver, A, b, x);  t_calc  = MPI_Wtime() - t_calc;

        /* Get info and release memory */
        HYPRE_Int num_iterations;
        HYPRE_Real final_res_norm;
        HYPRE_StructGMRESGetNumIterations( solver, &num_iterations );
        HYPRE_StructGMRESGetFinalRelativeResidualNorm( solver, &final_res_norm );

        if (my_pid == 0) {
            printf("\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Time cost %.5f %.5f %.5f %d\n", t_setup, t_calc, t_setup + t_calc, num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }

        // Destroy solver and preconditioner
        HYPRE_StructGMRESDestroy(solver);
        if      (strcmp(precond_name, "SMG" ) == 0) HYPRE_StructSMGDestroy (precond);
        else if (strcmp(precond_name, "PFMG") == 0) HYPRE_StructPFMGDestroy(precond);
    }

    // 计算真实残差
    HYPRE_StructVectorCopy(b, y);// y = b
    HYPRE_StructMatrixMatvec(-1.0, A, x, 1.0, y);// y += -A*x
    HYPRE_Real r_dot;
    r_dot = hypre_StructKrylovInnerProd(y, y);
    if (my_pid == 0) printf("\033[1;35mtrue ||r||/||b||= %20.16e\033[0m\n", sqrt(r_dot/b_dot));

    /* Free memory */
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(b); HYPRE_StructVectorDestroy(x); HYPRE_StructVectorDestroy(y);
    }// test loop

    /* Finalize HYPRE */
    HYPRE_Finalize();

    /* Finalize MPI */
    MPI_Finalize();
    return 0;
}