#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
using namespace std;

#define CORES_PER_NODE 120

#define SOCKET_PER_NODE 2
#define NUMAS_PER_SOCKET 2
#define CORES_PER_NUMA 30
#define NUMA_OFFSET 31 // 该程序将每个numa视作仅有31个核

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
    
    int procs_per_socket = procs_per_node / SOCKET_PER_NODE;
    assert(procs_per_node % SOCKET_PER_NODE == 0);// 每个socket上的进程数必须相同

    int procs_per_numa = procs_per_socket / NUMAS_PER_SOCKET;
    assert(procs_per_socket % NUMAS_PER_SOCKET == 0);// 每个numa上的进程数必须相同

    cout << num_nodes << " nodes each with " << procs_per_node << " procs" << endl;
    cout << "    each socket with " << procs_per_socket << " procs" << endl;
    cout << "    each numa   with " << procs_per_numa   << " procs" << endl;

    int threads_per_proc = (CORES_PER_NUMA * NUMAS_PER_SOCKET * SOCKET_PER_NODE) / procs_per_node;
    printf("    each proc has %d threads\n", threads_per_proc);

    vector<int> loc_core_ids_in_numa(procs_per_numa, -1);
    for (int p = 0; p < procs_per_numa; p++) {
        loc_core_ids_in_numa[p] = threads_per_proc * p;
    }

    string filename = "r1_N" + to_string(num_nodes) + "P" + to_string(num_procs);
    for (int g = 0; g < num_groups; g++)
        filename = filename + "_" + to_string(first_nids[g]) + "-" + to_string(first_nids[g] + node_cnts[g] - 1);
    FILE * fp = fopen(filename.c_str(), "w+");
    // print rankfile, 
    for (int N = 0; N < num_nodes; N++) {// by each node
        for (int n = 0; n < NUMAS_PER_SOCKET * SOCKET_PER_NODE; n++) {// by each numa
            int core_offset = n * NUMA_OFFSET;// 避开每个numa的前两个核
            for (int p = 0; p < procs_per_numa; p++) {
                int glb_pid = N * procs_per_node + n * procs_per_numa + p;
                int first_core_id = core_offset + loc_core_ids_in_numa[p];
                int last_core_id = first_core_id + threads_per_proc - 1;

                fprintf(fp, "rank %d=%s slot=%d-%d\n", 
                    glb_pid, node_names[N].c_str(), first_core_id, last_core_id);
            }
        }
    }
    fclose(fp);

    return 0;
}