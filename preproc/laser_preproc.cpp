#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
using namespace std;

typedef int idx_t;
typedef double data_t;
const int num_diag = 7;

int main(int argc, char* argv[])
{
    const char * mat_name = argv[1];
    const char * rhs_name = argv[2];
    const char * guess_name = argv[3];
    const idx_t m = atoi(argv[4]);
    const idx_t m3 = m*m*m;
    const idx_t offset[num_diag] = {-m*m, -m, -1, 0, 1, m, m*m};

    FILE * fp = fopen(mat_name, "r");

    idx_t ibeg, iend, jbeg, jend;// 闭区间
    fscanf(fp, "%d %d %d %d", &ibeg, &iend, &jbeg, &jend);
    assert(ibeg == jbeg && ibeg == 0);
    assert(iend == jend && iend == m3 - 1);

    data_t * diag_SOA = new data_t [m3 * num_diag];
    for (int p = 0; p < m3 * num_diag; p++)
        diag_SOA[p] = 0.0;
    
    idx_t r, c;
    data_t v;
    while (fscanf(fp, "%d %d %lf", &r, &c, &v) != EOF) {
        assert(r >= 0 && r < m3);
        assert(c >= 0 && c < m3);
        idx_t dist = c - r;
        for (idx_t d = 0; d < num_diag; d++) {
            if (offset[d] == dist) {
                diag_SOA[d * m3 + r] = v;
                break;
            }
        }
    }
    fclose(fp);
    for (int id = 0; id < 7; id++) {// 写出矩阵的各条对角线
        char buf[100];
        sprintf(buf, "array_a.%d", id);
        fp = fopen(buf, "wb+");
        int size = fwrite(diag_SOA + id * m3, sizeof(data_t), m3, fp);
        assert(size == m3);
        fclose(fp);
    }

    // 右端项
    data_t * rhs = new data_t [m3];
    fp = fopen(rhs_name, "r");
    fscanf(fp, "%d %d", &ibeg, &iend);
    assert(ibeg == jbeg && iend == jend);
    while (fscanf(fp, "%d %lf", &r, &v) != EOF) {
        rhs[r] = v;
    }
    fclose(fp);
    {// 写出右端项
        fp = fopen("array_b", "wb+");
        idx_t size = fwrite(rhs, sizeof(data_t), m3, fp);
        assert(size == m3);
        fclose(fp);
    }
    // 初始解
    fp = fopen(guess_name, "r");
    idx_t tmp[4];
    fscanf(fp, "%d %d %d %d", tmp, tmp+1, tmp+2, tmp+3);
    assert(tmp[0] == m && tmp[1] == m && tmp[2] == m);
    assert(tmp[3] == 1);
    r = 0;
    while (fscanf(fp, "%lf", &v) != EOF) {
        rhs[r] = v;
        r++;
    }
    assert(r == m3);
    fclose(fp);
    {// 写出初始解
        fp = fopen("array_x", "wb+");
        idx_t size = fwrite(rhs, sizeof(data_t), m3, fp);
        assert(size == m3);
        fclose(fp);
    }

    delete diag_SOA;
    delete rhs;
    return 0;
}
