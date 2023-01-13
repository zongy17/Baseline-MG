#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
using namespace std;

#define CORES_PER_NODE 124
#define NUM_NUMA 4
#define CORES_PER_NUMA 31

int main(int argc, char ** argv)
{
    string node_pfx = "96.10.130.";
    int argc_cnt = 1;

    const int num_groups = atoi(argv[argc_cnt++]);

    int first_nids[num_groups], node_cnts[num_groups];
    for (int g = 0; g < num_groups; g++) {
        first_nids[g] = atoi(argv[argc_cnt++]);
        node_cnts[g] = atoi(argv[argc_cnt++]);
    }

    vector<string> node_names;
    cout << "nodes' names" << endl;
    for (int g = 0;g < num_groups; g++) {
        for (int i = 0; i < node_cnts[g]; i++) {
            node_names.push_back(node_pfx + to_string(first_nids[g] + i));
            cout << node_names[node_names.size() - 1] << endl;
        }
    }
    const int num_nodes = node_names.size();
    
    int num_procs = atoi(argv[argc_cnt++]);
    int procs_per_node = num_procs / num_nodes;

    assert(num_procs % num_nodes == 0);// 每个节点上的进程数必须相同
    
    int threads_per_proc = CORES_PER_NODE / procs_per_node;
    if (threads_per_proc <= 1) {
        printf("Error: t/p <= 1!\n");
        exit(1);
    } else if (threads_per_proc <= 3) {
        printf("Warning: t/p = %d, may cause load imbalance\n", threads_per_proc);
    }

    vector<int> beg_core_ids(procs_per_node), cnt_cores(procs_per_node, threads_per_proc);
    int remain = CORES_PER_NODE - threads_per_proc * procs_per_node;
    if (remain > 0) {
        assert(procs_per_node % NUM_NUMA == 0);
        int procs_per_numa = procs_per_node / NUM_NUMA;

        int p = 0;
        while (remain > 0) {
            assert(remain % NUM_NUMA == 0);
            for (int n = 0; n < NUM_NUMA; n++) {
                cnt_cores[n * procs_per_numa + p] ++;
                remain --;   
            }
            p++;
        }
    }

    beg_core_ids[0] = 0;
    for (int p = 1; p < procs_per_node; p++)
        beg_core_ids[p] = beg_core_ids[p-1] + cnt_cores[p-1];

    cout << num_nodes << " nodes each with " << procs_per_node << " procs" << endl;
    printf("    each proc has %d threads\n", threads_per_proc);

    string filename = "hck_N" + to_string(num_nodes) + "P" + to_string(num_procs);
    for (int g = 0; g < num_groups; g++)
        filename = filename + "_" + to_string(first_nids[g]) + "-" + to_string(first_nids[g] + node_cnts[g] - 1);

    FILE * fp = fopen(filename.c_str(), "w+");
    // print rankfile, 
    for (int N = 0; N < num_nodes; N++) {// by each node
        for (int p = 0; p < procs_per_node; p++) {
            int glb_pid = N * procs_per_node + p;
            int first_core_id = beg_core_ids[p];
            int last_core_id = first_core_id + cnt_cores[p] - 1;

            fprintf(fp, "rank %d=%s slot=%d-%d\n", 
                glb_pid, node_names[N].c_str(), first_core_id, last_core_id);
        }
    }
    fclose(fp);

    return 0;
}