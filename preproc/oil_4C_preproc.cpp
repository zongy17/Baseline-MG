#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>
#include <assert.h>
using namespace std;

#define DOF 4

struct TUPLE
{
    int r, c;
    double v;
    TUPLE(int r, int c, double v): r(r), c(c), v(v) {}
    bool operator < (const TUPLE & b) const {
        if      (r < b.r) return true;
        else if (r > b.r) return false;
        else              return c < b.c; 
    }
};

int main(int argc, char ** argv)
{
    const std::string pathname = std::string(argv[1]);
    int timeid = atoi(argv[2]);
    
    FILE * fp = nullptr;
    int ret;

    double v[DOF*DOF];

    vector<int> row_ptr, col_idx;
    vector<double> vals;

    fp = fopen((pathname + "/testA."+std::to_string(timeid)).c_str(), "r");
    int br, bc, last_br = -1, last_bc = -1;// block row, block col

    row_ptr.push_back(0);
    vector<TUPLE> curr4row;
    while (fscanf(fp, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &br, &bc, v, v+1, v+2, v+3, v+4, v+5, v+6, v+7, v+8, v+9, v+10, v+11, v+12, v+13, v+14, v+15) != EOF) {
        if (br == last_br) {// 仍为原来那4行
            assert(bc > last_bc);
        }
        else {// 从此开始新的4行
            assert(br > last_br);
            if (last_br != -1) {
                // 将前4行的内容进行排序，并追补到正式数据末尾
                int cnt[DOF];
                for (int i = 0; i < DOF; i++) 
                    cnt[i] = 0;

                sort(curr4row.begin(), curr4row.end());
                for (int i = 0; i < curr4row.size(); i++) {
                    if (curr4row[i].v != 0.0) {// 只输出非零元？
                        int r = curr4row[i].r;
                        cnt[r - last_br * DOF] ++;
                        col_idx.push_back(curr4row[i].c);
                        vals.push_back(curr4row[i].v);
                    }
                }
                for (int i = 0; i < DOF; i++)
                    row_ptr.push_back(row_ptr[row_ptr.size() - 1] + cnt[i]);
            }
            // 然后清空
            curr4row.clear();
#ifdef DEBUG
            for (int i = 0; i < row_ptr.size(); i++)
                printf(" %d", row_ptr[i]);
            printf("\n");
            for (int i = 0; i < col_idx.size(); i++)
                printf(" %d", col_idx[i]);
            printf("\n");
            for (int i = 0; i < vals.size(); i++)
                printf(" %d", vals[i]);
            printf("\n");
#endif
        }

        for (int i = 0; i < DOF; i++)
        for (int j = 0; j < DOF; j++) {
            int r = br * DOF + i;
            int c = bc * DOF + j;
            curr4row.push_back(TUPLE(r, c, v[i*DOF + j]));
        }

        last_br = br;
        last_bc = bc;
    }
    // 跳出循环时，curr4row中应该不为空
    if (curr4row.size() > 0) {
        int cnt[DOF];
        for (int i = 0; i < DOF; i++) 
            cnt[i] = 0;

        sort(curr4row.begin(), curr4row.end());
        for (int i = 0; i < curr4row.size(); i++) {
            int r = curr4row[i].r;
            cnt[r - last_br * DOF] ++;
            col_idx.push_back(curr4row[i].c);
            vals.push_back(curr4row[i].v);
        }
        for (int i = 0; i < DOF; i++)
            row_ptr.push_back(row_ptr[row_ptr.size() - 1] + cnt[i]);

        // 然后清空
        curr4row.clear();
    }
    fclose(fp);

    assert(row_ptr[0] == 0);
    int tot_nnz = row_ptr[row_ptr.size() - 1];
    assert(tot_nnz == col_idx.size());
    assert(tot_nnz == vals.size());
    int nrows = row_ptr.size() - 1;
    printf(" nrows %d\n", nrows);

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
    while (fscanf(fp, "%lf %lf %lf %lf", v, v+1, v+2, v+3) != EOF) {
        for (int i = 0; i < DOF; i++)
            b.push_back(v[i]);
    }
    fclose(fp);
    fp = fopen("b.bin", "wb");
    ret = fwrite(b.data(), sizeof(double), b.size(), fp); 
    // assert(ret == nrows);
    fclose(fp);

    printf("b : %.6e %.6e ... %.6e %.6e\n", b[0], b[1], b[b.size()-2], b[b.size()-1]);

    return 0;
}