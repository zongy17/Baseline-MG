#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <assert.h>
using namespace std;

int main(int argc, char ** argv)
{
    const std::string pathname = std::string(argv[1]);
    int timeid = atoi(argv[2]);
    
    FILE * fp = nullptr;
    vector<int> row_ptr, col_idx;
    vector<double> vals;

    fp = fopen((pathname + "/testA."+std::to_string(timeid)).c_str(), "r");
    int r, c, last_row = -1, last_col = -1;
    double v;
    while (fscanf(fp, "%d %d %lf", &r, &c, &v) != EOF) {
        if (r == last_row) {// 仍为原来那行
            assert(c > last_col);
        }
        else {// 从此开始新一行
            assert(r > last_row);
            row_ptr.push_back(col_idx.size());
        }
        col_idx.push_back(c);
        vals.push_back(v);

        last_row = r;
        last_col = c;
    }
    row_ptr.push_back(col_idx.size());// 最后一行
    fclose(fp);

    assert(row_ptr[0] == 0);
    int tot_nnz = row_ptr[row_ptr.size() - 1];
    assert(tot_nnz == col_idx.size());
    assert(tot_nnz == vals.size());
    int nrows = row_ptr.size() - 1;
    printf(" nrows %d\n", nrows);

    int ret;
    // 输出矩阵的二进制
    fp = fopen("Ai.bin", "wb");
    ret = fwrite(row_ptr.data(), sizeof(int), row_ptr.size(), fp); assert(ret == row_ptr.size());
    fclose(fp);
    fp = fopen("Aj.bin", "wb");
    ret = fwrite(col_idx.data(), sizeof(int), col_idx.size(), fp); assert(ret == tot_nnz);
    fclose(fp);
    fp = fopen("Av.bin", "wb");
    ret = fwrite(vals.data(), sizeof(double), vals.size(), fp); assert(ret == tot_nnz);
    fclose(fp);

    vector<double> b;
    fp = fopen((pathname + "/testb."+std::to_string(timeid)).c_str(), "r");
    while (fscanf(fp, "%lf", &v) != EOF)
        b.push_back(v);
    fclose(fp);
    fp = fopen("b.bin", "wb");
    ret = fwrite(b.data(), sizeof(double), b.size(), fp); assert(ret == row_ptr.size() - 1);
    fclose(fp);

    vector<double> x0;
    fp = fopen((pathname + "/testx0."+std::to_string(timeid)).c_str(), "r");
    while (fscanf(fp, "%lf", &v) != EOF)
        x0.push_back(v);
    fclose(fp);
    fp = fopen("x0.bin", "wb");
    ret = fwrite(x0.data(), sizeof(double), x0.size(), fp); assert(ret == row_ptr.size() - 1);
    fclose(fp);

    return 0;
}