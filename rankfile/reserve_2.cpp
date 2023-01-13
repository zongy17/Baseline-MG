#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
using namespace std;

#define CORES_PER_NODE 120

#define SOCKET_PER_NODE 2
#define NUMAS_PER_SOCKET 2
#define CORES_PER_NUMA 30 // 该生成rankfile的程序避开了每个numa的前两个核心（0和1）
#define CORE_ID_OFFSET 2

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

    char filename[100];
    sprintf(filename, "N%dP%d\0", num_nodes, num_procs);
    FILE * fp = fopen(filename, "w+");
    // print rankfile, 
    for (int N = 0; N < num_nodes; N++) {// by each node
        for (int n = 0; n < NUMAS_PER_SOCKET * SOCKET_PER_NODE; n++) {// by each numa
            int core_offset = (n+1) * CORE_ID_OFFSET;// 避开每个numa的前两个核
            for (int p = 0; p < procs_per_numa; p++) {
                int glb_pid = N * procs_per_node + n * procs_per_numa + p;
                int first_core_id = CORES_PER_NUMA * n + loc_core_ids_in_numa[p] + core_offset;
                fprintf(fp, "rank %d=%s slot=%d-%d\n", 
                    glb_pid, node_names[N].c_str(), first_core_id, first_core_id + threads_per_proc - 1);
            }
        }
    }
    fclose(fp);

    return 0;
}