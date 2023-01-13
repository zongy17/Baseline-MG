#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <assert.h>
using namespace std;

#ifndef USE_BigInt
typedef int idx_t;
#else
typedef long long idx_t;
#endif

typedef float data_t;

const idx_t num_diag = 19;
const idx_t jik_off[num_diag * 3] = {
    -1, -1,  0,// 0
    -1,  0, -1,// 1
    -1,  0,  0,// 2
    -1,  0,  1,// 3
    -1,  1,  0,// 4

    0, -1, -1,// 5
    0, -1,  0,// 6
    0, -1,  1,// 7
    0,  0, -1,// 8
    0,  0,  0,// 9
    0,  0,  1,// 10
    0,  1, -1,// 11
    0,  1,  0,// 12
    0,  1,  1,// 13

    1, -1,  0,// 14
    1,  0, -1,
    1,  0,  0,
    1,  0,  1,
    1,  1,  0 // 18
};

double dot(data_t * x, data_t * y, idx_t n) {
    double tmp = 0.0;
    for (idx_t i = 0; i < n; i++)
        tmp += (double) x[i] * (double) y[i];
    return tmp;
}

const idx_t map_storage_to_memory_diag[19] = {
 // 0  1   2  3   4   5  6  7  8   9 10  11 12  13  14 15  16 17  18 
    9, 6, 12, 2, 16, 18, 4, 0, 14, 8, 5, 11, 1, 15, 10, 7, 13, 3, 17
};

int main(int argc, char ** argv)
{
    const std::string pathname = std::string(argv[1]);
    const int timeid = atoi(argv[2]);

    idx_t lon, lat, lev;
    if (pathname.find("50km") != string::npos) {
        lon =  720; lat =  360;
    } else if (pathname.find("25km") != string::npos) {
        lon = 1440; lat =  720;
    } else if (pathname.find("12_5km") != string::npos) {
        lon = 2880; lat = 1440;
        if (sizeof(idx_t) != 8) {
            printf("Error: 12.5km must use Long Long instead of int\n");
            exit(1);
        }
    } else exit(1);
    lev = 62;

    idx_t tot_elem = lon * lat * lev;
    idx_t tot_nnz = tot_elem * num_diag;
    data_t * vec_buf = new data_t [tot_elem];
    data_t * mat_buf = new data_t [tot_nnz];
    for (idx_t i = 0; i < tot_nnz; i++) {
        mat_buf[i] = i;
    }

    FILE * fp = nullptr; idx_t ret;

    for (idx_t d = 0; d < num_diag; d++) {// 逐条对角线读入A
        printf("reading from %s ", (pathname + "/array_a.new" + to_string(d+1)).c_str());
        fp = fopen((pathname + "/array_a.new" + to_string(d+1)).c_str(), "rb"); assert(fp);
        ret = fread(vec_buf, sizeof(data_t), tot_elem, fp); assert(ret == tot_elem);
        fclose(fp);
        printf(" done\n");

        const idx_t target = map_storage_to_memory_diag[d];
        for (idx_t j = 0; j < lat; j++)
        for (idx_t i = 0; i < lon; i++)
        for (idx_t k = 0; k < lev; k++) {
            idx_t src_id = (j * lev + k) * lon + i;
            idx_t trg_id = (j * lon + i) * lev + k;
            assert(src_id < tot_elem && trg_id < tot_elem);
            mat_buf[trg_id * num_diag + target] = vec_buf[src_id];
        }
    }

    vector<idx_t> row_ptr, col_idx;
    vector<data_t> vals;
    assert(lon % 2 == 0);
    idx_t half_lon = lon / 2;
    row_ptr.push_back(0);
    for (idx_t j = 0; j < lat; j++)
    for (idx_t i = 0; i < lon; i++)
    for (idx_t k = 0; k < lev; k++) {
        idx_t cnt = 0;
        data_t * ptr = mat_buf + ((j * lon + i) * lev + k) * num_diag;
        for (idx_t d = 0; d < num_diag; d++) {
            idx_t ngb_j = j + jik_off[3*d  ];
            idx_t ngb_i = i + jik_off[3*d+1];
            idx_t ngb_k = k + jik_off[3*d+2];
            if (ngb_k < 0 || ngb_k >= lev) {
                assert(ptr[d] == 0.0);
                continue;
            }
            cnt ++;
            // 首先将可能超出[0, lon)范围内的ngb_i限制回来
            if      (ngb_i <  0  ) ngb_i += lon;
            else if (ngb_i >= lon) ngb_i -= lon;
            // 然后使用跨极点规律
            if (ngb_j < 0 || ngb_j >= lat) {
                ngb_j = j;
                if (ngb_i < half_lon) ngb_i += half_lon;
                else                  ngb_i -= half_lon;
            }
            idx_t col = (ngb_j * lon + ngb_i) * lev + ngb_k;
            assert(col >= 0 && col < tot_elem);
            col_idx.push_back(col);
            vals.push_back(ptr[d]);
        }
        idx_t last_nnz = row_ptr[row_ptr.size() - 1];
        row_ptr.push_back(last_nnz + cnt);
    }

    assert(row_ptr[0] == 0);
    idx_t real_nnz = row_ptr[row_ptr.size() - 1]; assert(real_nnz == idx_t(col_idx.size())); assert(real_nnz == idx_t(vals.size()));
    idx_t nrows = row_ptr.size() - 1; assert(nrows == tot_elem);
    printf(" nrows %lld\n", nrows);

    data_t * x = new data_t [tot_elem];
    {// 读入x
        fp = fopen((pathname + "/array_x.new" + to_string(timeid)).c_str(), "rb"); assert(fp);
        ret = fread(vec_buf, sizeof(data_t), tot_elem, fp); assert(ret == tot_elem);
        fclose(fp);

        for (idx_t j = 0; j < lat; j++)
        for (idx_t i = 0; i < lon; i++)
        for (idx_t k = 0; k < lev; k++) {
            idx_t src_id = (j * lev + k) * lon + i;
            idx_t trg_id = (j * lon + i) * lev + k;
            assert(src_id < tot_elem && trg_id < tot_elem);
            x[trg_id] = vec_buf[src_id];    
        }
    }

    data_t * b = new data_t [tot_elem];
    {// 读入b
        fp = fopen((pathname + "/array_b.new" + to_string(timeid)).c_str(), "rb"); assert(fp);
        ret = fread(vec_buf, sizeof(data_t), tot_elem, fp); assert(ret == tot_elem);
        fclose(fp);

        for (idx_t j = 0; j < lat; j++)
        for (idx_t i = 0; i < lon; i++)
        for (idx_t k = 0; k < lev; k++) {
            idx_t src_id = (j * lev + k) * lon + i;
            idx_t trg_id = (j * lon + i) * lev + k;
            assert(src_id < tot_elem && trg_id < tot_elem);
            b[trg_id] = vec_buf[src_id];
        }
    }
    // 检验：向量的点积
    double thedot = 0.0;
    thedot = dot(x, x, tot_elem);
    printf("( x,  x) = %.15e\n", thedot);
    thedot = dot(b, b, tot_elem);
    printf("( b,  b) = %.15e\n", thedot);

    // 输出向量的二进制
    fp = fopen("b.bin", "wb");
    ret = fwrite(b, sizeof(data_t), tot_elem, fp); assert(ret == nrows);
    fclose(fp);
    fp = fopen("x0.bin", "wb");
    ret = fwrite(x, sizeof(data_t), tot_elem, fp); assert(ret == nrows);
    fclose(fp);

    // 检验：A*x，要在输出b之后再做
    for (idx_t i = 0; i < nrows; i++) {
        idx_t pbeg = row_ptr[i], pend = row_ptr[i+1];
        data_t tmp = 0.0;
        for (idx_t p = pbeg; p < pend; p++)
            tmp += vals[p] * x[col_idx[p]];
        b[i] = tmp;
    }
    thedot = dot(b, b, tot_elem);
    printf("(Ax, Ax) = %.15e\n", thedot);

    // 输出矩阵的二进制
    fp = fopen("Ai.bin", "wb");
    ret = fwrite(row_ptr.data(), sizeof(idx_t), row_ptr.size(), fp); assert(ret == idx_t(row_ptr.size()));
    fclose(fp);
    fp = fopen("Aj.bin", "wb");
    ret = fwrite(col_idx.data(), sizeof(idx_t), col_idx.size(), fp); assert(ret == real_nnz);
    fclose(fp);
    fp = fopen("Av.bin", "wb");
    ret = fwrite(vals.data(), sizeof(data_t), vals.size(), fp); assert(ret == real_nnz);
    fclose(fp);

    delete vec_buf;
    delete mat_buf;
    return 0;
}