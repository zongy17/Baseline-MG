#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
using namespace std;

#define CORES_PER_NODE 128

int main(int argc, char ** argv)
{
    string first_node_name = string(argv[1]);
    int first_node_idx = atoi(argv[2]);
    int num_nodes = atoi(argv[3]);

    vector<string> node_names;
    cout << "nodes' names" << endl;
    for (int i = 0; i < num_nodes; i++) {
        node_names.push_back(first_node_name + to_string(first_node_idx + i));
        cout << node_names[node_names.size() - 1] << endl;
    }
    
    int num_procs = atoi(argv[4]);
    int need_spare = atoi(argv[5]);
    if (need_spare) printf("spare first %d cores each numa.\n", need_spare);

    int procs_per_node = num_procs / num_nodes;
    assert(num_procs % num_nodes == 0);// 每个节点上的进程数必须相同
    
    assert(CORES_PER_NODE % procs_per_node == 0);
    int threads_per_proc = CORES_PER_NODE / procs_per_node;

    cout << num_nodes << " nodes each with " << procs_per_node << " procs" << endl;
    printf("    each proc has %d threads\n", threads_per_proc);

    vector<int> beg_core_ids(procs_per_node);
    for (int p = 0; p < procs_per_node; p++) {
        beg_core_ids[p] = threads_per_proc * p;
    }

    char filename[100];
    sprintf(filename, "N%dP%d\0", num_nodes, num_procs);
    FILE * fp = fopen(filename, "w+");
    // print rankfile, 
    for (int N = 0; N < num_nodes; N++) {// by each node
        for (int p = 0; p < procs_per_node; p++) {
            int glb_pid = N * procs_per_node + p;
            int first_core_id = beg_core_ids[p];
            int last_core_id = first_core_id + threads_per_proc - 1;
            if (need_spare && first_core_id % 32 == 0) first_core_id += need_spare;
            fprintf(fp, "rank %d=%s slot=%d-%d\n", 
                glb_pid, node_names[N].c_str(), first_core_id, last_core_id);
        }
    }
    fclose(fp);

    return 0;
}