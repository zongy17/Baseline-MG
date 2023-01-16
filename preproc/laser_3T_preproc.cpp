#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
using namespace std;

typedef int idx_t;
typedef double data_t;
const int num_diag = 7;
const int num_dof = 3;

int main(int argc, char* argv[])
{
    const char * mat_name = argv[1];
    const char * rhs_name = argv[2];
    const idx_t nx = atoi(argv[3]);
    const idx_t ny = atoi(argv[4]);
    const idx_t nz = atoi(argv[5]);
    const idx_t tot_elems = nx * ny * nz;
    const idx_t offset[num_diag] = {-nz*nx, -nz, -1, 0, 1, nz, nz*nx};
    printf("nx %d ny %d nz %d\n", nx, ny, nz);

    FILE * fp = fopen(mat_name, "r");

    idx_t ibeg, iend, jbeg, jend;// 闭区间
    fscanf(fp, "%d %d %d %d", &ibeg, &iend, &jbeg, &jend);
    assert(ibeg == jbeg && ibeg == 0);
    assert(iend == jend && iend == tot_elems*num_dof - 1);

    data_t * diag_SOA = new data_t [tot_elems * num_diag * num_dof * num_dof];
    for (int p = 0; p < tot_elems * num_diag * num_dof * num_dof; p++)
        diag_SOA[p] = 0.0;
    
    idx_t r, c;
    data_t v;
    while (fscanf(fp, "%d %d %lf", &r, &c, &v) != EOF) {
        assert(r >= 0 && r < tot_elems * num_dof);
        assert(c >= 0 && c < tot_elems * num_dof);
        idx_t blk_r = r % tot_elems, var_r = r / tot_elems;
        idx_t blk_c = c % tot_elems, var_c = c / tot_elems;
        idx_t dist = blk_c - blk_r;

        for (idx_t d = 0; d < num_diag; d++) {
            if (offset[d] == dist) {
                data_t * ptr = diag_SOA + (d * tot_elems + blk_r) * num_dof*num_dof;
                ptr[var_r * num_dof + var_c] = v;
                break;
            }
        }
    }
    fclose(fp);
    for (int id = 0; id < num_diag; id++) {// 写出矩阵的各条对角线
        char buf[100];
        sprintf(buf, "array_a.%d", id);
        fp = fopen(buf, "wb+");
        int size = fwrite(diag_SOA + id * tot_elems * num_dof * num_dof, sizeof(data_t), tot_elems * num_dof * num_dof, fp);
        assert(size == tot_elems * num_dof * num_dof);
        fclose(fp);
    }

    // 右端项
    data_t * rhs = new data_t [tot_elems * num_dof];
    fp = fopen(rhs_name, "r");
    fscanf(fp, "%d %d", &ibeg, &iend);
    assert(ibeg == jbeg && iend == jend);
    while (fscanf(fp, "%d %lf", &r, &v) != EOF) {
        idx_t blk = r % tot_elems;
        idx_t var = r / tot_elems;
        rhs[blk * num_dof + var] = v;
    }
    double sum = 0.0;
    for (int i = 0; i < tot_elems * num_dof; i++)
        sum += rhs[i] * rhs[i];
    printf("(b, b) = %.10e\n", sum);

    fclose(fp);
    {// 写出右端项
        fp = fopen("array_b", "wb+");
        idx_t size = fwrite(rhs, sizeof(data_t), tot_elems * num_dof, fp);
        assert(size == tot_elems * num_dof);
        fclose(fp);
    }

    delete diag_SOA;
    delete rhs;
    return 0;
}
