/*
Dat Nguyen
331007653
Compile: g++ -o algo_project.exe ./algo_project.cc
Run: algo_project.exe
*/

#include <iostream>
#include <fstream>
#include <random>
#include <list>
#include <assert.h>
#include <algorithm>
#include <tuple>
#include <cstring>
#include <unordered_map>
#include <queue>
#include <sys/time.h>
#include <bitset>

using namespace std;

#define GRAPH_SIZE 5000
#define TEST_NUM 5
#define ST_PAIRS 5

class Graph {
private:
    int _size;
    list<tuple<int, int>>**_adj_lists;

    void delete_tuple_in_list(list<tuple<int, int>>* L, int v) {
        for (list<tuple<int, int>>::iterator it = L->begin(); it != L->end(); ++it) {
            if (get<0>(*it) == v) {
                L->erase(it);
                break;
            }
        }
    }
public:
    Graph(int size): _size(size), _adj_lists((list<tuple<int, int>>**)malloc(size * sizeof (list<tuple<int, int>>*))) {
        for (int i = 0; i < _size; i++) {
            _adj_lists[i] = new list<tuple<int, int>>();
        }
    }

    ~Graph() {
        for (int i = 0; i < _size; i++) {
            delete _adj_lists[i];
        }
        free(_adj_lists);
    }

    int size() { return _size; }

    list<tuple<int, int>>* adj_list(int v) { return _adj_lists[v]; }

    void add_edge(int u, int v ,int weight) {
        _adj_lists[u]->push_back(tuple<int, int>(v, weight));
        _adj_lists[v]->push_back(tuple<int, int>(u, weight));
    }

    void delete_edge(int u, int v) {
        // _adj_lists[u]->remove(v);
        // _adj_lists[v]->remove(u);
        delete_tuple_in_list(_adj_lists[u], v);
        delete_tuple_in_list(_adj_lists[v], u);
    }

    bool is_adj(int u, int v) {
        list<tuple<int, int>>* adj_list = _adj_lists[u];
        for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
            if (get<0>(*it) == v) {
                return true;
            }
        }
        return false;
    }

    int weight(int u, int v) {
        // printf("u %d, v %d\n", u, v);
        list<tuple<int, int>>* adj_list = _adj_lists[u];
        assert(adj_list != NULL);
        for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
            if (get<0>(*it) == v) {
                return get<1>(*it);
            }
        }
        return -1;
    }

    void print_graph(int print_size=-1) {
        print_size = print_size == -1 ? _size : min(print_size, _size);
        for (int i = 0; i < print_size; i++) {
            list<tuple<int, int>>* adj_list = _adj_lists[i];
            printf("[%d]", i);
            for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
                printf("\t(%d, %d)", get<0>(*it), get<1>(*it));
            }
            printf("\n");
        }
    }
};

class MaxHeap {
private:
    int _size;
    int _max_size;
    int* _vertices; // store the vertices aka the heap
    int* _values; // map each vertice to its value
    int* _positions;    // map each vertice to its pos

    int _dad_pos(int pos)         { return (pos - 1) / 2; }
    int left_child_pos(int pos)     { return (2 * pos) + 1; }
    int right_child_pos(int pos)    { return (2 * pos) + 2; }

    int _dad(int v)               { return _vertices[_dad_pos(_positions[v])]; }
    int left_child(int v)           { return _vertices[left_child_pos(_positions[v])]; }
    int right_child(int v)          { return _vertices[right_child_pos(_positions[v])]; }

    bool pos_is_leaf(int pos) {
        return (pos > (_size / 2) && pos <= _size);
    }

    bool is_leaf(int v) {
        return pos_is_leaf(_positions[v]);
    }

    bool is_root(int v) {
        return _positions[v] == 0;
    }

    void pos_swap (int pos1, int pos2) {
        assert(pos1 < _size);
        assert(pos2 < _size);
        int tmp;
        tmp = _vertices[pos1];
        _vertices[pos1] = _vertices[pos2];
        _vertices[pos2] = tmp;
    }

    void swap (int u, int v) {
        if (u == v) return;
        // u and v are vertices name
        int u_pos = _positions[u];
        int v_pos = _positions[v];
        assert(u_pos < _size && u_pos != -1);
        assert(v_pos < _size && v_pos != -1);

        int tmp;
        tmp = _vertices[u_pos];
        assert(tmp == u);
        _vertices[u_pos] = _vertices[v_pos];
        _vertices[v_pos] = tmp;

        assert(_vertices[u_pos] == v);
        assert(_vertices[v_pos] == u);

        _positions[u] = v_pos;
        _positions[v] = u_pos;

    }

    void heapify(int subtree) {
        if (is_leaf(subtree) || _positions[subtree] >= _size || _positions[subtree] < 0) return;
        int left = left_child(subtree);
        int right = right_child(subtree);
        if (_values[subtree] < _values[left] || _values[subtree] < _values[right]) {
 
            if (_values[left] > _values[right]) {
                swap(subtree, left);
                heapify(left);
            }
            else {
                swap(subtree, right);
                heapify(right);
            }
        }
    }

public:
    MaxHeap(int max_size=5000) :
        _size(0),
        _max_size(max_size),
        _vertices((int*)malloc(max_size * sizeof(int))),
        _values((int*)malloc(max_size * sizeof(int))),
        _positions((int*)malloc(max_size * sizeof(int)))
    {
        for (int i = 0; i < max_size; i++) {
            _vertices[i]    = -1;
            _values[i]      = -1;
            _positions[i]   = -1;
        }
    }

    int max_size() { return _max_size; }

    bool is_in(int v) { return _positions[v] != -1; }

    bool is_empty() { return _vertices[0] == -1;}

    ~MaxHeap() {
        free(_vertices);
        free(_values);
        free(_positions);
    }

    int MAXIMUM() {
        // return vertice with max value, might return -1
        return _vertices[0];
    }

    void INSERT(int v, int value) {
        if (_positions[v] != -1) {
            printf("Already in heap \n");
            return;
        }
        assert(v < _max_size);
        _positions[v] = _size;
        _vertices[_size] = v;
        _values[v] = value;
        // int current = _size;
        _size++;
        while(!is_root(v) && _values[v] > _values[_dad(v)]) {
            swap(v, _dad(v));
        }
    }

    void DELETE(int v) {
        //delete element at pos 
        assert(v < _max_size);
        assert(_positions[v] != -1);
        int pos_v = _positions[v];
        // _vertices[pos] = _vertices[--_size];
        
        int last = _vertices[_size-1];
        assert(_positions[last] != -1);
        swap(v, last);
        _positions[v] = -1;
        _values[v] = -1;
        _vertices[_size-1] = -1;
        _size -= 1;
        if (v != last) {
            assert(_vertices[pos_v] == last);
        }
        heapify(last);
    }

    void print(int limit = -1) {
        if (limit == -1) {
            limit = _max_size;
        } else {
            limit = min(limit, _max_size);
        }
        assert(limit >= 0);
        printf("Vertices: [");
        for (int i = 0; i < limit; i++) {
            printf(" %d", _vertices[i]);
        }
        printf("]\n");

        printf("Values: [");
        for (int i = 0; i < limit; i++) {
            printf(" %d", _values[i]);
        }
        printf("]\n");

        printf("Pos: [");
        for (int i = 0; i < limit; i++) {
            printf(" %d", _positions[i]);
        }
        printf("]\n");
    }
};

class DisjointSet {
    int* _dad;
    int* _rank;
    int _max_size;
 
public:
    DisjointSet(int max_size) :
            _max_size(max_size),
            _dad((int*)malloc(max_size * sizeof(int))),
            _rank((int*)malloc(max_size * sizeof(int))){
        for (int i = 0; i < max_size; i++) {
            _dad[i] = -1;
            _rank[i] = -1;
        }
    }

    ~DisjointSet() {
        free(_dad);
        free(_rank);
    }

    void MakeSet(int v)
    {
        _dad[v] = v;
        _rank[v] = 0;
    }
 
    int Find(int v)
    {   
        if (_dad[v] == -1) return _dad[v];
        if (_dad[v] != v)
        {
            // path compression
            _dad[v] = Find(_dad[v]);
        }
 
        return _dad[v];
    }
 
    // Perform Union of two subsets
    void Union(int a, int b) {
        int x = Find(a);
        int y = Find(b);
 
        if (x == y) {
            return;
        }
 
        if (_rank[x] > _rank[y]) {
            _dad[y] = x;
        }
        else if (_rank[x] < _rank[y]) {
            _dad[x] = y;
        }
        else {
            _dad[x] = y;
            _rank[y]++;
        }
    }

    bool of_same_set(int a, int b) {
        return Find(a) == Find(b);
    }

    int find_or_add(int v) {
        // find operation, if does not exist, add v as new set
        int find = Find(v);
        if (find == -1) {
            MakeSet(v);
            return v;
        }
        return find;
    }
};

class VEMapper {
private:
    int _bits;
    int _bits_mask;
    int _maximum;

public: 
    VEMapper(int largest_v = 4999) : _bits(0), _maximum(0) {
        while (largest_v > 0) {
            _bits += 1;
            largest_v = largest_v >> 1;
        }
        _bits_mask = (1 << _bits) -1;
        // printf("_bits %d\n", _bits);
        // printf("_bits_mask %d\n", _bits_mask);
        // std::cout<<std::bitset<32>(_bits_mask)<<std::endl;
        int smaller = 0;
        int larger = 1;
        int e = encode(smaller, larger);
        int new_smaller = 0;
        int new_larger = 0;
        decode(e, new_smaller, new_larger);
        assert(smaller == new_smaller);
        assert(larger == new_larger);
        // printf("_bits", _bits);
        // printf("_bits_mask %d", _bits_mask);
    }

    int bits() { return _bits; }

    int maximum()  { return _maximum; }

    int encode(int u, int v) {
        assert(u != v);
        assert(u < _bits_mask);
        assert(v < _bits_mask);
        int small = min(u, v);
        int large = max(u, v);
        int diff = large-small;
        diff = diff << _bits;
        diff = diff | small;
        // std::cout<<std::bitset<32>(diff)<<std::endl;
        return diff;
    }

    void decode(int e, int& smaller, int& larger) {
        // std::cout<<std::bitset<32>(e)<<std::endl;
        int diff = e >> _bits;
        smaller = e & _bits_mask;
        // std::cout<<std::bitset<32>(smaller)<<std::endl;
        larger = smaller + diff;
        // std::cout<<std::bitset<32>(larger)<<std::endl;
    }

    void preview(int e) {
        _maximum = max(_maximum, e);
    }
};

Graph* create_g1(int random_seed=42) {
    printf("Creating G1\n");
    int graph_size = GRAPH_SIZE;
    int average_vertex_deg = 6;
    int target_num_edges = average_vertex_deg * graph_size;
    int num_edges_count = 0;

    Graph* G = new Graph(graph_size /*graph size*/);
    
    int weight = 0;
    // make a cycle of vertices first
    for (int v = 0; v < graph_size; v++) {
        int w = (v + 1) % graph_size;
        weight = (rand() + 1) % (INT32_MAX -1);
        G->add_edge(v, w, weight);
        num_edges_count += 1;
    }

    // add the rest of edges randomly until reach target
    srand(random_seed);
    while (num_edges_count != target_num_edges) {
        int v, w;
        do {
            v = rand() % graph_size;
            do {
                w = rand() % graph_size;
            } while (w == v);

        } while (G->is_adj(v, w));
        
        weight = (rand() + 1) % (INT32_MAX -1);
        G->add_edge(v, w, weight);
        num_edges_count += 1;
    }
    printf("G1 created\n");
    printf("==============================\n");
    return G;
}

Graph* create_g2(int random_seed=42) {
    printf("Creating G2\n");
    int graph_size = GRAPH_SIZE;
    int adj_percent = 20; // each v is adj to roughly 20% of other nodes
    int target_adj = (graph_size * 20 / 100) - 1; 

    Graph* G = new Graph(graph_size /*graph size*/);

    int weight = 0;
    // make a cycle of vertices first
    for (int v = 0; v < graph_size; v++) {
        int w = (v + 1) % graph_size;
        weight = (rand() + 1) % (INT32_MAX -1);
        G->add_edge(v, w, weight);
    }

    // for each vertex, add edges until target adj is reached
    srand(random_seed);
    int progress_step = 100;
    for (int v = 0; v < graph_size; v++) {
        int adj_count = G->adj_list(v)->size();
        while (adj_count < target_adj) {
            int w;
            do {
                w = rand() % graph_size;
            } while (w == v || G->is_adj(v, w));
            weight = (rand() + 1) % (INT32_MAX -1);
            G->add_edge(v, w, weight);
            adj_count += 1;
        }

        if (v % progress_step == 0) {
            int perc = v * 100 / graph_size;
            printf("[%.*s", perc, "==================================================================================================");
            printf(">");
            printf("%.*s]", 100 - perc, "                                                                                                     ");
            printf(" %d%%\n", perc);
        }
    }

    printf("G2 created\n");
    printf("==============================\n");
    return G;
}

// =====================================================================

enum STATUS {
    UNSEEN,
    INTREE,
    FRINGER
};

void MaxBanwidthPath_0(Graph* G, int s, int t, int* bw, int* dad) {
    // Dijkstra without heap structure

    int graph_size = G->size();
    assert(s < graph_size);
    assert(t < graph_size);


    STATUS status  [graph_size];
    for (int v = 0; v < graph_size; v++) {
        status[v] = UNSEEN;
        bw[v] = 0;
        dad[v] = -1;
    }
    status[s] = INTREE; bw[s] = -1 /* bw = -1 = inf */; dad[s] = -1;

    list<tuple<int, int>>* adj_list = NULL;
    adj_list = G->adj_list(s);
    for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
        int w = get<0>(*it);
        status[w] = FRINGER;
        bw[w] = get<1>(*it);
        dad[w] = s;
    }

    int fringer_w_max_bw = -1;
    do {
        fringer_w_max_bw = -1;
        int max_bw = -1;
        for (int i = 0; i < graph_size; i ++) {
            if (status[i] == FRINGER && bw[i] > max_bw) {
                fringer_w_max_bw = i;
                max_bw = bw[i];
            }
        }

        if (fringer_w_max_bw != -1) {
            int v = fringer_w_max_bw;
            status[v] = INTREE;
            adj_list = G->adj_list(v);
            for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
                int w = get<0>(*it);
                if (status[w] == UNSEEN) {
                    status[w] = FRINGER;
                    dad[w] = v;
                    bw[w] = min(bw[v], get<1>(*it));
                } else if (status[w] == FRINGER && bw[w] < min(bw[v], get<1>(*it))) {
                    dad[w] = v;
                    bw[w] = min(bw[v], get<1>(*it));
                }
            }
        }

    } while (fringer_w_max_bw != -1);
}

void MaxBanwidthPath_1(Graph* G, int s, int t, int* bw, int* dad) {
    // Dijkstra with heap structure

    int graph_size = G->size();
    assert(s < graph_size);
    assert(t < graph_size);


    STATUS status  [graph_size];
    for (int v = 0; v < graph_size; v++) {
        status[v] = UNSEEN;
        bw[v] = 0;
        dad[v] = -1;
    }
    status[s] = INTREE; bw[s] = -1 /* bw = -1 = inf */; dad[s] = -1;

    MaxHeap heap(G->size());

    list<tuple<int, int>>* adj_list = NULL;
    adj_list = G->adj_list(s);
    for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
        int w = get<0>(*it);
        status[w] = FRINGER;
        bw[w] = get<1>(*it);
        dad[w] = s;
        heap.INSERT(w, bw[w]);
    }

    int fringer_w_max_bw = -1;
    do {
        fringer_w_max_bw = -1;

        fringer_w_max_bw = heap.MAXIMUM();

        if (fringer_w_max_bw != -1) {
            int v = fringer_w_max_bw;
            heap.DELETE(v);
            status[v] = INTREE;
            adj_list = G->adj_list(v);
            for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
                int w = get<0>(*it);
                if (status[w] == UNSEEN) {
                    status[w] = FRINGER;
                    dad[w] = v;
                    bw[w] = min(bw[v], get<1>(*it));
                    heap.INSERT(w, bw[w]);
                } else if (status[w] == FRINGER && bw[w] < min(bw[v], get<1>(*it))) {
                    heap.DELETE(w);
                    dad[w] = v;
                    bw[w] = min(bw[v], get<1>(*it));
                    heap.INSERT(w, bw[w]);
                }
            }
        }

    } while (fringer_w_max_bw != -1);
}

void MaxBanwidthPath_2(Graph* G, int s, int t, int* bw, int* dad) {
    // Kruskal with edges sorted by heap sort
    int graph_size = G->size();
    assert(s < graph_size);
    assert(t < graph_size);

    // note, change edge heap size base of vemapper bits
    VEMapper* ve_mapper = new VEMapper();
    // printf("VE mapper created\n");

    list<tuple<int, int>>* adj_list = NULL;
    // loop over all edges and find maximum edge idx

    int progress_step = 100;
    for (int i = 0; i < graph_size; i++) {
        int u = i;
        adj_list = G->adj_list(u);
        for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
            int v = get<0>(*it);
            int e = ve_mapper->encode(u, v);
            ve_mapper->preview(e);
        }
        // if (u % progress_step == 0) {
        //     int perc = u * 100 / graph_size;
        //     printf("[%.*s", perc, "==================================================================================================");
        //     printf(">");
        //     printf("%.*s]", 100 - perc, "                                                                                                     ");
        //     printf(" %d%%\n", perc);
        // }
    }
    // printf("MaxHeap max size = %d\n", ve_mapper->maximum());

    MaxHeap* edge_heap = new MaxHeap(ve_mapper->maximum());


    // loop over all edges and add them to heap
    // printf("Adding all edges to max heap!\n");
    for (int i = 0; i < graph_size; i++) {
        int u = i;
        adj_list = G->adj_list(u);
        for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
            int v = get<0>(*it);
            int weight = get<1>(*it);
            int e = ve_mapper->encode(u, v);
            if (!edge_heap->is_in(e)) edge_heap->INSERT(e, weight);
        }
        // if (u % progress_step == 0) {
        //     int perc = u * 100 / graph_size;
        //     printf("[%.*s", perc, "==================================================================================================");
        //     printf(">");
        //     printf("%.*s]", 100 - perc, "                                                                                                     ");
        //     printf(" %d%%\n", perc);
        // }
    }


    // printf("Adding edges to tree structure!\n");
    Graph* tree = new Graph(graph_size);
    DisjointSet* disjoint_set = new DisjointSet(graph_size);
    while (edge_heap->MAXIMUM() != -1) {
        int e = edge_heap->MAXIMUM();
        int u = 0; int v = 0;
        ve_mapper->decode(e, u, v);
        int u_set = disjoint_set->find_or_add(u);
        int v_set = disjoint_set->find_or_add(v);
        if (u_set != v_set) {
            tree->add_edge(u, v, 1);
            disjoint_set->Union(u, v);
        }
        edge_heap->DELETE(e);
    }

    // // // perform BFS
    // printf("Perform BFS\n");
    queue<int> Q;
    int visited [GRAPH_SIZE] = {};
    Q.push(s);
    while (!Q.empty()) {
        int v = Q.front();
        Q.pop();
        adj_list = tree->adj_list(v);
        visited[v] = 1;
        for (list<tuple<int, int>>::iterator it = adj_list->begin(); it != adj_list->end(); ++it) {
            int w = get<0>(*it);
            // int weight = get<1>(*it);
            if (visited[w] == 0) {
                dad[w] = v;
                bw[w] = G->weight(w, v);
                if (w == t) goto RET;
                Q.push(w);
            }
        }
    }

    RET:
    delete tree;
    delete ve_mapper;
    delete edge_heap;
    delete disjoint_set;
}

// =====================================================================
void print_arr(string prefix, int* arr, int size) {
    printf("%s [", prefix.c_str());
    for(int i = 0; i < size; i ++) {
        printf(" %d", arr[i]);
        
    }
    printf(" ]\n");
}

void my_test() {
    Graph* G = new Graph(10 /*graph size*/);
    G->add_edge(1, 2, 3);
    G->add_edge(4, 5, 6);
    // printf("1, 2 is adj %d\n", G->is_adj(1, 2));
    // printf("1, 3 is adj %d\n", G->is_adj(1, 3));
    G->delete_edge(1, 2);
    G->print_graph(10);
    G->weight(2, 1);
    G->weight(8, 9);

    // MaxHeap* heap = new MaxHeap((1<<27)-1);
    // heap->INSERT(1, 12);
    // heap->INSERT(2, 6);
    // heap->INSERT(3, 6);
    // heap->INSERT(2, 8);
    // // heap->print();

    // heap->INSERT(4, 18);
    // // heap->print();
    // heap->DELETE(1);
    // // heap->print();

    // VEMapper* ve_mapper = new VEMapper();
    // int e = ve_mapper->encode(12, 3000);
    // printf("Encode: %d\n", e);

}

void TEST() {
    srand(time(NULL));

    double runtimes0 [2][TEST_NUM] = {};
    double runtimes1 [2][TEST_NUM] = {};
    double runtimes2 [2][TEST_NUM] = {};

    int bw[GRAPH_SIZE];
    int dad[GRAPH_SIZE];

    ofstream outfile;
    outfile.open("outfile.txt");
    outfile << "Output: \n";

    for (int i = 0; i < TEST_NUM; i++) {
        printf("Test %d\n", i);
        double total_elapsed_G1;
        Graph* G1 = create_g1(time(NULL));
        Graph* G2 = create_g2(time(NULL));
        for (int j = 0; j < ST_PAIRS; j++) {
            int s = rand() % GRAPH_SIZE;
            int t = rand() % GRAPH_SIZE;

            // test algo 0
            printf("Algo 0\n");
            {
                // G1
                printf("G1\n");
                memset(bw, 0, GRAPH_SIZE * sizeof(int));
                memset(dad, 0, GRAPH_SIZE * sizeof(int));
                struct timeval begin, end;
                gettimeofday(&begin, 0);
                MaxBanwidthPath_0(G1, s, t, bw, dad);
                gettimeofday(&end, 0);
                double elapsed = (end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec)*1e-6;
                runtimes0[0][i] += elapsed;
            }
            {
                // G2
                printf("G2\n");
                memset(bw, 0, GRAPH_SIZE * sizeof(int));
                memset(dad, 0, GRAPH_SIZE * sizeof(int));
                struct timeval begin, end;
                gettimeofday(&begin, 0);
                MaxBanwidthPath_0(G2, s, t, bw, dad);
                gettimeofday(&end, 0);
                double elapsed = (end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec)*1e-6;
                runtimes0[1][i] += elapsed;
            }

            // test algo 1
            printf("Algo 1\n");
            {
                // G1
                printf("G1\n");
                memset(bw, 0, GRAPH_SIZE * sizeof(int));
                memset(dad, 0, GRAPH_SIZE * sizeof(int));
                struct timeval begin, end;
                gettimeofday(&begin, 0);
                MaxBanwidthPath_1(G1, s, t, bw, dad);
                gettimeofday(&end, 0);
                double elapsed = (end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec)*1e-6;
                runtimes1[0][i] += elapsed;
            }
            {
                // G2
                printf("G2\n");
                memset(bw, 0, GRAPH_SIZE * sizeof(int));
                memset(dad, 0, GRAPH_SIZE * sizeof(int));
                struct timeval begin, end;
                gettimeofday(&begin, 0);
                MaxBanwidthPath_1(G2, s, t, bw, dad);
                gettimeofday(&end, 0);
                double elapsed = (end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec)*1e-6;
                runtimes1[1][i] += elapsed;
            }

            // test algo 2
            printf("Algo 2\n");
            {
                // G1
                printf("G1\n");
                memset(bw, 0, GRAPH_SIZE * sizeof(int));
                memset(dad, 0, GRAPH_SIZE * sizeof(int));
                struct timeval begin, end;
                gettimeofday(&begin, 0);
                MaxBanwidthPath_2(G1, s, t, bw, dad);
                gettimeofday(&end, 0);
                double elapsed = (end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec)*1e-6;
                runtimes2[0][i] += elapsed;
            }
            {
                // G2
                printf("G2\n");
                memset(bw, 0, GRAPH_SIZE * sizeof(int));
                memset(dad, 0, GRAPH_SIZE * sizeof(int));
                struct timeval begin, end;
                gettimeofday(&begin, 0);
                MaxBanwidthPath_2(G2, s, t, bw, dad);
                gettimeofday(&end, 0);
                double elapsed = (end.tv_sec - begin.tv_sec) + (end.tv_usec - begin.tv_usec)*1e-6;
                runtimes2[1][i] += elapsed;
            }
        }


        delete G1;
        delete G2;
    }
    outfile.close();
}

int main () {
    TEST();
    // my_test();
    // Graph* G1 = create_g1();
    // int bw[5000];
    // int dad[5000];
    // MaxBanwidthPath_0(G1, 1, 2000, bw, dad);
    // MaxBanwidthPath_1(G1, 1, 2000, bw, dad);
    // MaxBanwidthPath_2(G1, 1, 2000, bw, dad);
    // // G1->print_graph(10);
    // // print_arr("bandwidth", bw, 5000);
    // // print_arr("dad", dad, 5000);
    // Graph* G2 = create_g2();
    // MaxBanwidthPath_0(G2, 1, 2000, bw, dad);
    // MaxBanwidthPath_1(G2, 1, 2000, bw, dad);

}