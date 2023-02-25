#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "HYPRE_sstruct_ls.h"
#include "_hypre_sstruct_ls.h"

#define MIN(a, b) (a) < (b) ? (a) : (b)
#ifndef TEST_CNT
#define TEST_CNT 1
#endif

int main(int argc, char * argv[])
{
    int my_pid, num_procs;
    HYPRE_Int glb_nrows;

    MPI_Init(&argc, &argv);

    int cnt = 1;
    const char * case_name = argv[cnt++];
    HYPRE_Int gx = atoi(argv[cnt++]),
        gy = atoi(argv[cnt++]),
        gz = atoi(argv[cnt++]);
    HYPRE_Int px = atoi(argv[cnt++]),
        py = atoi(argv[cnt++]),
        pz = atoi(argv[cnt++]);
    const char * precond_name = argv[cnt++];
    if (strcmp(case_name, "LASER-3T") == 0) {
        glb_nrows = gx * gy * gz * 3;
    }
    else if (strcmp(case_name, "SOLID") == 0) {
        glb_nrows = gx * gy * gz * 3;
    }
    else {
        MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
        if (my_pid == 0) printf("Error: Unrecognized case name %s\n", case_name);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    HYPRE_Int lx, ly, lz;
    if (gx % px != 0 || gy % py != 0 || gz % pz != 0) {// check if input is legal
        MPI_Comm_rank(MPI_COMM_WORLD, &my_pid);
        if (my_pid == 0) printf("Error: Illegal input! gx/y/z must be fully divided by px/y/z!\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    } else {
        lx = gx / px;
        ly = gy / py;
        lz = gz / pz;
    }
    const HYPRE_Int nparts = 1;
    const HYPRE_Int nvars = 3;
    const HYPRE_Int part = 0;// part id
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
    HYPRE_SStructGrid     grid;
    HYPRE_SStructGraph    graph;
    HYPRE_SStructStencil  stencils[3];
    HYPRE_SStructMatrix   A;
    HYPRE_SStructVector   b, x, y;
    HYPRE_SStructSolver   solver, precond;

    for (HYPRE_Int test = 0; test < TEST_CNT; test++) {
    {// Setup a grid
        HYPRE_SStructGridCreate(cart_comm, ndim, nparts, &grid);// Create an empty 3D grid object
        HYPRE_SStructGridSetExtents(grid, part, ilower, iupper);// Add a new box to the grid
        HYPRE_SStructVariable var_types[nvars] = {
            HYPRE_SSTRUCT_VARIABLE_CELL, HYPRE_SSTRUCT_VARIABLE_CELL, HYPRE_SSTRUCT_VARIABLE_CELL
        };
        // Set the variable type and number of variables on each part.
        for (HYPRE_Int i = 0; i < nparts; i++)
            HYPRE_SStructGridSetVariables(grid, i, nvars, var_types);
        HYPRE_SStructGridAssemble(grid);// a collective call finalizing the grid assembly.
    }

    if (strcmp(case_name, "LASER-3T") == 0) {
        const HYPRE_Int num_diag = 7;
        const HYPRE_Int stencil_size = num_diag + 2;// 7 for 所有邻居结构点（含自己）的该变量，2 for 自己结构点的另两个变量
        HYPRE_Int offsets[stencil_size][ndim] = {
            {0,0,-1}, {0,-1,0}, {-1,0,0}, {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1},
            {0,0,0}, {0,0,0}
        };
        {// Define the discretization stencil
            for (HYPRE_Int s = 0; s < nvars; s++) {
                HYPRE_SStructStencilCreate(ndim, stencil_size, &stencils[s]);
                HYPRE_Int vars[stencil_size];
                // 前七个
                for (HYPRE_Int j = 0; j < num_diag; j++)
                    vars[j] = s;
                // 后两个
                if      (s == 0) { vars[num_diag] = 1; vars[num_diag+1] = 2; }
                else if (s == 1) { vars[num_diag] = 0; vars[num_diag+1] = 2; }
                else if (s == 2) { vars[num_diag] = 0; vars[num_diag+1] = 1; }
                else MPI_Abort(cart_comm, -20230101);
                for (HYPRE_Int e = 0; e < stencil_size; e++)
                    HYPRE_SStructStencilSetEntry(stencils[s], e, offsets[e], vars[e]);
            }
        }

        HYPRE_Int object_type;
        {// Setup the graph
            HYPRE_SStructGraphCreate(cart_comm, grid, &graph);
            object_type = HYPRE_SSTRUCT;
            HYPRE_SStructGraphSetObjectType(graph, object_type);
            // Assign the stencils[s] we created to variable[s]
            for (HYPRE_Int s = 0; s < nvars; s++)
                HYPRE_SStructGraphSetStencil(graph, part, s, stencils[s]);
            // Assemble the graph
            HYPRE_SStructGraphAssemble(graph);
        }

        HYPRE_Int stencil_indices[stencil_size];
        for (HYPRE_Int j = 0; j < stencil_size; j++)
            stencil_indices[j] = j;

        const HYPRE_Int nnz_len = 3 * (num_diag-1) + 3*3; assert(nnz_len == nvars * (num_diag+2));
        {// Set up a Struct Matrix and Struct Vector
            const HYPRE_Int num_elems = lx * ly * lz;
            HYPRE_Int tot_len = nnz_len * num_elems;
            HYPRE_Real * buf = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * tot_len);
            // 首先MPI-IO读入矩阵数据
            MPI_File fh = MPI_FILE_NULL;// 文件句柄
            MPI_Datatype etype = (sizeof(HYPRE_Real) == 4) ? MPI_FLOAT : MPI_DOUBLE;
            MPI_Datatype read_type = MPI_DATATYPE_NULL;// 写出类型
            HYPRE_Int size[4], subsize[4], start[4];
            size   [0] = gy       ; size   [1] = gx       ; size   [2] = gz       ; size   [3] = nnz_len;
            subsize[0] = ly       ; subsize[1] = lx       ; subsize[2] = lz       ; subsize[3] = nnz_len;
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
            HYPRE_SStructMatrixCreate(cart_comm, graph, &A);// Create an empty matrix object
            HYPRE_SStructMatrixSetObjectType(A, object_type);
            HYPRE_SStructMatrixInitialize(A);// Indicate that the matrix coefficients are ready to be set 
            HYPRE_Real * mat_0_0 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * num_diag * num_elems);
            HYPRE_Real * mat_1_1 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * num_diag * num_elems);
            HYPRE_Real * mat_2_2 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * num_diag * num_elems);
            HYPRE_Real * mat_0_1 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) *            num_elems);
            HYPRE_Real * mat_0_2 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) *            num_elems);
            HYPRE_Real * mat_1_0 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) *            num_elems);
            HYPRE_Real * mat_1_2 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) *            num_elems);
            HYPRE_Real * mat_2_0 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) *            num_elems);
            HYPRE_Real * mat_2_1 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) *            num_elems);
            {// 格式转换
                for (HYPRE_Int i = 0; i < num_elems; i++) {
                    const HYPRE_Real * src_ptr = buf + nnz_len * i;
                    {// 第一个变量
                        HYPRE_Real * dst_ptr = mat_0_0 + i * num_diag;
                        // 跟别的邻居结构点的第一个变量
                        dst_ptr[0] = src_ptr[0];  dst_ptr[1] = src_ptr[3];  dst_ptr[2] = src_ptr[6];
                        dst_ptr[3] = src_ptr[9];  mat_0_1[i] = src_ptr[10]; mat_0_2[i] = src_ptr[11];
                        dst_ptr[4] = src_ptr[18]; dst_ptr[5] = src_ptr[21]; dst_ptr[6] = src_ptr[24];
                    }
                    {// 第二个变量
                        HYPRE_Real * dst_ptr = mat_1_1 + i * num_diag;
                        // 跟别的邻居结构点的第二个变量
                        dst_ptr[0]= src_ptr[1];  dst_ptr[1]= src_ptr[4];   dst_ptr[2]= src_ptr[7];
                        mat_1_0[i]= src_ptr[12]; dst_ptr[3]= src_ptr[13];  mat_1_2[i]= src_ptr[14];
                        dst_ptr[4]= src_ptr[19]; dst_ptr[5]= src_ptr[22];  dst_ptr[6]= src_ptr[25];
                    }
                    {// 第三个变量
                        HYPRE_Real * dst_ptr = mat_2_2 + i * num_diag;
                        // 跟别的邻居结构点的第三个变量
                        dst_ptr[0]= src_ptr[2];  dst_ptr[1]= src_ptr[5];  dst_ptr[2]= src_ptr[8];
                        mat_2_0[i]= src_ptr[15]; mat_2_1[i]= src_ptr[16]; dst_ptr[3]= src_ptr[17]; 
                        dst_ptr[4]= src_ptr[20]; dst_ptr[5]= src_ptr[23]; dst_ptr[6]= src_ptr[26];  
                    }
                }
            }

            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 0, num_diag, stencil_indices             , mat_0_0);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 0,   1     , stencil_indices + num_diag  , mat_0_1);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 0,   1     , stencil_indices + num_diag+1, mat_0_2);

            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 1, num_diag, stencil_indices             , mat_1_1);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 1,   1     , stencil_indices + num_diag  , mat_1_0);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 1,   1     , stencil_indices + num_diag+1, mat_1_2);

            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 2, num_diag, stencil_indices             , mat_2_2);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 2,   1     , stencil_indices + num_diag  , mat_2_0);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 2,   1     , stencil_indices + num_diag+1, mat_2_1);

            HYPRE_SStructMatrixAssemble(A);// a collective call finalizing the matrix assembly.

            MPI_Type_free(&read_type);
            free(buf);
            free(mat_0_0); free(mat_0_1); free(mat_0_2);
            free(mat_1_1); free(mat_1_0); free(mat_1_2);
            free(mat_2_2); free(mat_2_0); free(mat_2_1);

            // 向量
            size   [3] = 3;// #dof
            subsize[3] = 3;
            start  [3] = 0;
            MPI_Type_create_subarray(4, size, subsize, start, MPI_ORDER_C, etype, &read_type);
            MPI_Type_commit(&read_type);

            HYPRE_SStructVectorCreate(cart_comm, grid, &b);// Create an empty vector object
            HYPRE_SStructVectorCreate(cart_comm, grid, &x);
            HYPRE_SStructVectorCreate(cart_comm, grid, &y);
            HYPRE_SStructVectorSetObjectType(b, object_type);// Set the object type for the vectors to be the same as was already set for the matrix
            HYPRE_SStructVectorSetObjectType(x, object_type);
            HYPRE_SStructVectorSetObjectType(y, object_type);
            HYPRE_SStructVectorInitialize(b);// Indicate that the vector coefficients are ready to be set
            HYPRE_SStructVectorInitialize(x);
            HYPRE_SStructVectorInitialize(y);

            tot_len = 3 * num_elems;
            buf = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * tot_len);
            HYPRE_Real * vec_0 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * num_elems);
            HYPRE_Real * vec_1 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * num_elems);
            HYPRE_Real * vec_2 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * num_elems);
            {// b
                ret = MPI_File_open(cart_comm, "b.AOS.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
                if (ret != MPI_SUCCESS) {
                    printf("Could not open file: ret %d\n", ret);
                }
                MPI_File_set_view(fh, displacement, etype, read_type, "native", MPI_INFO_NULL);
                MPI_File_read_all(fh, buf, tot_len, etype, &status);
                MPI_File_close(&fh);
                for (HYPRE_Int i = 0; i < num_elems; i++) {
                    vec_0[i] = buf[i*3  ];
                    vec_1[i] = buf[i*3+1];
                    vec_2[i] = buf[i*3+2];
                }
                HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, 0, vec_0);
                HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, 1, vec_1);
                HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, 2, vec_2);
            }
            {// x
                ret = MPI_File_open(cart_comm, "x.AOS.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
                if (ret != MPI_SUCCESS) {
                    printf("Could not open file: ret %d\n", ret);
                }
                MPI_File_set_view(fh, displacement, etype, read_type, "native", MPI_INFO_NULL);
                MPI_File_read_all(fh, buf, tot_len, etype, &status);
                MPI_File_close(&fh);
                for (HYPRE_Int i = 0; i < num_elems; i++) {
                    vec_0[i] = buf[i*3  ];
                    vec_1[i] = buf[i*3+1];
                    vec_2[i] = buf[i*3+2];
                }
                HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, 0, vec_0);
                HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, 1, vec_1);
                HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, 2, vec_2);
            }
            HYPRE_SStructVectorAssemble(b);
            HYPRE_SStructVectorAssemble(x);
            HYPRE_SStructVectorAssemble(y);// 对y也要创建
            
            MPI_Type_free(&read_type);
            free(buf);
            free(vec_0); free(vec_1); free(vec_2);
        }
    }
    else if (strcmp(case_name, "SOLID") == 0) {
        const HYPRE_Int read_num_diag = 15, pfmg_num_diag = 27;
        const HYPRE_Int stencil_size = pfmg_num_diag * nvars;
        HYPRE_Int translate_3d15_to_3d27[read_num_diag] = {
            0, 1, 3, 4, 9, 10, 12, 13, 14, 16, 17, 22, 23, 25, 26
        };
        HYPRE_Int offsets[stencil_size][ndim] = {
            // (k, i, j)
            // 对第一个变量
                                        {-1,-1,-1},{0,-1,-1},                           {1,-1,-1},
                                        {-1,0,-1},{0,0,-1},                             {1, 0,-1},
                                                                    {-1,1,-1},{ 0,1,-1},{1,1,-1},
                                        {-1,-1,0},{0,-1,0},                             {1,-1,0},
                                        {-1,0,0}, {0,0,0},{1,0,0},
            {-1,1,0},                             {0,1,0},{1,1,0},
            {-1,-1,1},{0,-1,1},{1,-1,1},
            {-1,0,1},                             {0,0,1},{1,0,1},
            {-1,1,1},                             {0,1,1},{1,1,1},
            // 对第二个变量
                                        {-1,-1,-1},{0,-1,-1},                           {1,-1,-1},
                                        {-1,0,-1},{0,0,-1},                             {1, 0,-1},
                                                                    {-1,1,-1},{ 0,1,-1},{1,1,-1},
                                        {-1,-1,0},{0,-1,0},                             {1,-1,0},
                                        {-1,0,0}, {0,0,0},{1,0,0},
            {-1,1,0},                             {0,1,0},{1,1,0},
            {-1,-1,1},{0,-1,1},{1,-1,1},
            {-1,0,1},                             {0,0,1},{1,0,1},
            {-1,1,1},                             {0,1,1},{1,1,1},
            // 对第三个变量
                                        {-1,-1,-1},{0,-1,-1},                           {1,-1,-1},
                                        {-1,0,-1},{0,0,-1},                             {1, 0,-1},
                                                                    {-1,1,-1},{ 0,1,-1},{1,1,-1},
                                        {-1,-1,0},{0,-1,0},                             {1,-1,0},
                                        {-1,0,0}, {0,0,0},{1,0,0},
            {-1,1,0},                             {0,1,0},{1,1,0},
            {-1,-1,1},{0,-1,1},{1,-1,1},
            {-1,0,1},                             {0,0,1},{1,0,1},
            {-1,1,1},                             {0,1,1},{1,1,1}
        };
        {// Define the discretization stencil
            for (HYPRE_Int s = 0; s < nvars; s++) {
                HYPRE_SStructStencilCreate(ndim, stencil_size, &stencils[s]);
                for (HYPRE_Int e = 0; e < stencil_size; e++) {
                    HYPRE_Int ngb_var = e / pfmg_num_diag;
                    HYPRE_SStructStencilSetEntry(stencils[s], e, offsets[e], ngb_var);
                }
            }
        }

        HYPRE_Int object_type;
        {// Setup the graph
            HYPRE_SStructGraphCreate(cart_comm, grid, &graph);
            object_type = HYPRE_SSTRUCT;
            HYPRE_SStructGraphSetObjectType(graph, object_type);
            // Assign the stencils[s] we created to variable[s]
            for (HYPRE_Int s = 0; s < nvars; s++)
                HYPRE_SStructGraphSetStencil(graph, part, s, stencils[s]);
            // Assemble the graph
            HYPRE_SStructGraphAssemble(graph);
        }

        HYPRE_Int stencil_indices[stencil_size];
        for (HYPRE_Int j = 0; j < stencil_size; j++)
            stencil_indices[j] = j;

        const HYPRE_Int lms = nvars * nvars;
        {// Set up a Struct Matrix and Struct Vector
            const HYPRE_Int num_elems = lx * ly * lz;
            HYPRE_Int tot_len = num_elems * read_num_diag * lms;
            HYPRE_Real * buf = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * tot_len);
            // 首先MPI-IO读入矩阵数据
            MPI_File fh = MPI_FILE_NULL;// 文件句柄
            MPI_Datatype etype = (sizeof(HYPRE_Real) == 4) ? MPI_FLOAT : MPI_DOUBLE;
            MPI_Datatype read_type = MPI_DATATYPE_NULL;// 写出类型
            HYPRE_Int size[4], subsize[4], start[4];
            size   [0] = gy       ; size   [1] = gx       ; size   [2] = gz       ; size   [3] = read_num_diag * lms;
            subsize[0] = ly       ; subsize[1] = lx       ; subsize[2] = lz       ; subsize[3] = read_num_diag * lms;
            start  [0] = pidy * ly; start  [1] = pidx * lx; start  [2] = pidz * lz; start  [3] = 0                  ;

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
            HYPRE_SStructMatrixCreate(cart_comm, graph, &A);// Create an empty matrix object
            HYPRE_SStructMatrixSetObjectType(A, object_type);
            HYPRE_SStructMatrixInitialize(A);// Indicate that the matrix coefficients are ready to be set 
            const size_t bytes = sizeof(HYPRE_Real) * pfmg_num_diag * num_elems;
            HYPRE_Real * mat_0_0 = (HYPRE_Real *) malloc (bytes);
            HYPRE_Real * mat_0_1 = (HYPRE_Real *) malloc (bytes);
            HYPRE_Real * mat_0_2 = (HYPRE_Real *) malloc (bytes);
            HYPRE_Real * mat_1_0 = (HYPRE_Real *) malloc (bytes);
            HYPRE_Real * mat_1_1 = (HYPRE_Real *) malloc (bytes);
            HYPRE_Real * mat_1_2 = (HYPRE_Real *) malloc (bytes);
            HYPRE_Real * mat_2_0 = (HYPRE_Real *) malloc (bytes);
            HYPRE_Real * mat_2_1 = (HYPRE_Real *) malloc (bytes);
            HYPRE_Real * mat_2_2 = (HYPRE_Real *) malloc (bytes);
            for (size_t i = 0; i < pfmg_num_diag * num_elems; i++) {
                mat_0_0[i] = 0.0; mat_0_1[i] = 0.0; mat_0_2[i] = 0.0;
                mat_1_0[i] = 0.0; mat_1_1[i] = 0.0; mat_1_2[i] = 0.0;
                mat_2_0[i] = 0.0; mat_2_1[i] = 0.0; mat_2_2[i] = 0.0;
            }
            {// 格式转换
                assert(lms == 9);
                for (HYPRE_Int i = 0; i < num_elems; i++) {
                    for (HYPRE_Int d = 0; d < read_num_diag; d++) {// 每个对角块有9个元素
                        const HYPRE_Real * src_ptr = buf + (i * read_num_diag + d) * lms;
                        const HYPRE_Int loc = i * pfmg_num_diag + translate_3d15_to_3d27[d];
                        mat_0_0[loc] = src_ptr[0]; mat_0_1[loc] = src_ptr[1]; mat_0_2[loc] = src_ptr[2];
                        mat_1_0[loc] = src_ptr[3]; mat_1_1[loc] = src_ptr[4]; mat_1_2[loc] = src_ptr[5];
                        mat_2_0[loc] = src_ptr[6]; mat_2_1[loc] = src_ptr[7]; mat_2_2[loc] = src_ptr[8];
                    }
                }
            }

            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 0, pfmg_num_diag, stencil_indices                  , mat_0_0);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 0, pfmg_num_diag, stencil_indices + pfmg_num_diag  , mat_0_1);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 0, pfmg_num_diag, stencil_indices + pfmg_num_diag*2, mat_0_2);

            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 1, pfmg_num_diag, stencil_indices                  , mat_1_0);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 1, pfmg_num_diag, stencil_indices + pfmg_num_diag  , mat_1_1);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 1, pfmg_num_diag, stencil_indices + pfmg_num_diag*2, mat_1_2);

            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 2, pfmg_num_diag, stencil_indices                  , mat_2_0);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 2, pfmg_num_diag, stencil_indices + pfmg_num_diag  , mat_2_1);
            HYPRE_SStructMatrixSetBoxValues(A, part, ilower, iupper, 2, pfmg_num_diag, stencil_indices + pfmg_num_diag*2, mat_2_2);

            HYPRE_SStructMatrixAssemble(A);// a collective call finalizing the matrix assembly.

            MPI_Type_free(&read_type);
            free(buf);
            free(mat_0_0); free(mat_0_1); free(mat_0_2);
            free(mat_1_1); free(mat_1_0); free(mat_1_2);
            free(mat_2_2); free(mat_2_0); free(mat_2_1);

            // 向量
            size   [3] = nvars;// #dof
            subsize[3] = nvars;
            start  [3] = 0;
            MPI_Type_create_subarray(4, size, subsize, start, MPI_ORDER_C, etype, &read_type);
            MPI_Type_commit(&read_type);

            HYPRE_SStructVectorCreate(cart_comm, grid, &b);// Create an empty vector object
            HYPRE_SStructVectorCreate(cart_comm, grid, &x);
            HYPRE_SStructVectorCreate(cart_comm, grid, &y);
            HYPRE_SStructVectorSetObjectType(b, object_type);// Set the object type for the vectors to be the same as was already set for the matrix
            HYPRE_SStructVectorSetObjectType(x, object_type);
            HYPRE_SStructVectorSetObjectType(y, object_type);
            HYPRE_SStructVectorInitialize(b);// Indicate that the vector coefficients are ready to be set
            HYPRE_SStructVectorInitialize(x);
            HYPRE_SStructVectorInitialize(y);

            tot_len = nvars * num_elems;
            buf = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * tot_len);
            HYPRE_Real * vec_0 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * num_elems);
            HYPRE_Real * vec_1 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * num_elems);
            HYPRE_Real * vec_2 = (HYPRE_Real *) malloc (sizeof(HYPRE_Real) * num_elems);
            {// b
                ret = MPI_File_open(cart_comm, "b.AOS.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
                if (ret != MPI_SUCCESS) {
                    printf("Could not open file: ret %d\n", ret);
                }
                MPI_File_set_view(fh, displacement, etype, read_type, "native", MPI_INFO_NULL);
                MPI_File_read_all(fh, buf, tot_len, etype, &status);
                MPI_File_close(&fh);
                for (HYPRE_Int i = 0; i < num_elems; i++) {
                    vec_0[i] = buf[i*3  ];
                    vec_1[i] = buf[i*3+1];
                    vec_2[i] = buf[i*3+2];
                }
                HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, 0, vec_0);
                HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, 1, vec_1);
                HYPRE_SStructVectorSetBoxValues(b, part, ilower, iupper, 2, vec_2);
            }
            {// x
                ret = MPI_File_open(cart_comm, "x.AOS.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
                if (ret != MPI_SUCCESS) {
                    printf("Could not open file: ret %d\n", ret);
                }
                MPI_File_set_view(fh, displacement, etype, read_type, "native", MPI_INFO_NULL);
                MPI_File_read_all(fh, buf, tot_len, etype, &status);
                MPI_File_close(&fh);
                for (HYPRE_Int i = 0; i < num_elems; i++) {
                    vec_0[i] = buf[i*3  ];
                    vec_1[i] = buf[i*3+1];
                    vec_2[i] = buf[i*3+2];
                }
                HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, 0, vec_0);
                HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, 1, vec_1);
                HYPRE_SStructVectorSetBoxValues(x, part, ilower, iupper, 2, vec_2);
            }
            
            HYPRE_SStructVectorAssemble(b);
            HYPRE_SStructVectorAssemble(x);
            HYPRE_SStructVectorAssemble(y);// 对y也要创建
            
            MPI_Type_free(&read_type);
            free(buf);
            free(vec_0); free(vec_1); free(vec_2);
        }
    }

    HYPRE_Real b_dot;
    b_dot = hypre_SStructKrylovInnerProd(b, b);
    if (my_pid == 0) printf("(  b,   b) = %.20e\n",  b_dot);
    {// Check for DEBUG purpose
        HYPRE_Real x_dot, Ab_dot;
        x_dot = hypre_SStructKrylovInnerProd(x, x);
        if (my_pid == 0) printf("(  x,   x) = %.20e\n",  x_dot);
        HYPRE_SStructMatrixMatvec(1.0, A, b, 0.0, y);
        Ab_dot= hypre_SStructKrylovInnerProd(y, y);
        if (my_pid == 0) printf("(A*b, A*b) = %.20e\n", Ab_dot);
    }

    
    {// Setup a struct solver, preconditioner and RUN
        HYPRE_SStructPCGCreate(cart_comm, &solver);
        HYPRE_SStructPCGSetMaxIter(solver, 1000);
        HYPRE_SStructPCGSetTol(solver, 1.0e-09);

        HYPRE_SStructPCGSetTwoNorm(solver, 1 );
        HYPRE_SStructPCGSetRelChange(solver, 0 );
        HYPRE_SStructPCGSetPrintLevel(solver, 2);
        HYPRE_SStructPCGSetLogging(solver, 1);

        if (strcmp(precond_name, "SysPFMG") == 0) {
            if (my_pid == 0) printf("  using SysPFMG\n");
            HYPRE_SStructSysPFMGCreate(cart_comm, &precond);
            HYPRE_SStructSysPFMGSetTol(precond, 0.0);
            HYPRE_SStructSysPFMGSetMaxIter(precond, 1);
            HYPRE_SStructSysPFMGSetPrintLevel(precond, 2);

            HYPRE_SStructSysPFMGSetRelaxType(precond, 2);// 0: 不收敛！ 1: 284次 2: 59次
            HYPRE_SStructSysPFMGSetNumPreRelax(precond, 1);
            HYPRE_SStructSysPFMGSetNumPostRelax(precond, 1);
            HYPRE_SStructSysPFMGSetZeroGuess(precond);
            // HYPRE_SStructSysPFMGSetSkipRelax(precond, 1);
            
            HYPRE_SStructPCGSetPrecond(solver,  HYPRE_SStructSysPFMGSolve,
                                                HYPRE_SStructSysPFMGSetup, precond);
        }
        else if (strcmp(precond_name, "NOPC") == 0) { 
            if (my_pid == 0) printf("  NO precond!\n");
        }
        else {
            if (my_pid == 0) printf("  Error: Uncognized PC name: %s!\n", precond_name);
            MPI_Abort(MPI_COMM_WORLD, -101);
        }

        HYPRE_Real t_setup = MPI_Wtime();
        HYPRE_SStructPCGSetup(solver, A, b, x ); t_setup = MPI_Wtime() - t_setup;
        HYPRE_Real t_calc  = MPI_Wtime();
        HYPRE_SStructPCGSolve(solver, A, b, x);  t_calc  = MPI_Wtime() - t_calc;

        // Get info and release memory
        HYPRE_Int num_iterations;
        HYPRE_Real final_res_norm;
        HYPRE_SStructPCGGetNumIterations( solver, &num_iterations );
        HYPRE_SStructPCGGetFinalRelativeResidualNorm( solver, &final_res_norm );

        if (my_pid == 0) {
            printf("\n");
            printf("Iterations = %d\n", num_iterations);
            printf("Time cost %.5f %.5f %.5f %d\n", t_setup, t_calc, t_setup + t_calc, num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }

        // Destroy solver and preconditioner
        HYPRE_SStructPCGDestroy(solver);
        if (strcmp(precond_name, "sysPFMG") == 0) HYPRE_SStructSysPFMGDestroy(precond);
    }

    // 计算真实残差
    HYPRE_SStructVectorCopy(b, y);// y = b
    HYPRE_SStructMatrixMatvec(-1.0, A, x, 1.0, y);// y += -A*x
    HYPRE_Real r_dot;
    r_dot = hypre_SStructKrylovInnerProd(y, y);
    if (my_pid == 0) printf("\033[1;35mtrue ||r||/||b||= %20.16e\033[0m\n", sqrt(r_dot/b_dot));
    

    /* Free memory */
    HYPRE_SStructGridDestroy(grid);
    HYPRE_SStructGraphDestroy(graph);
    for (HYPRE_Int s = 0; s < nvars; s++)
        HYPRE_SStructStencilDestroy(stencils[s]);
    HYPRE_SStructMatrixDestroy(A);
    HYPRE_SStructVectorDestroy(b); HYPRE_SStructVectorDestroy(x); HYPRE_SStructVectorDestroy(y);
    }

    /* Finalize HYPRE */
    HYPRE_Finalize();

    /* Finalize MPI */
    MPI_Finalize();
    return 0;
}