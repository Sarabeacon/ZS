

#ifndef bclist_parallel_h
#define bclist_parallel_h

#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <sstream>
#include <utility>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <time.h>
#include <cmath>
#include <omp.h>



#define FREESUB
#define COUNT_ONLY

//#define IS_DEBUGGING

#define TRIM_BY_TWO_HOP_NEIGHBORS

#define USE_CORE_REDUCTION


//#define TRIM_BY_TWO_HOP_NEIGHBORS

#define SZ(x) ((int)x.size())
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

#define NTWOHOPS 10000

#define SAMPLE_RATE 0.01

#define MAX_N 10000000
#define MAX_K 20


using namespace std;

typedef unsigned long long lint;



class VertexDegree{
public:
    unsigned vertex;
    unsigned degree;
    inline bool operator < (const VertexDegree &other) const {
        return degree > other.degree || (degree == other.degree && vertex < other.vertex);
    }
    
    VertexDegree();
    VertexDegree(unsigned v, unsigned d);
    ~VertexDegree();
};

inline bool cmpByDegAsc(VertexDegree &v1, VertexDegree &v2){
    return v1.degree < v2.degree || (v1.degree == v2.degree && v1.vertex < v2.vertex);
}

/*typedef struct VertexRank {
    unsigned vertex;
    unsigned rank;
} VertexRank;*/

class Tools{
public:
    static vector<unsigned> intersection(vector<unsigned> &vec1, vector<unsigned> &vec2);
    static vector<unsigned> intersection(vector<unsigned> &vec1, vector<unsigned> &vec2, int offset1, int offset2);
    static unsigned intersection_count(vector<unsigned> &vec1, vector<unsigned> &vec2);
    static lint choose(lint n, lint k);
    static lint my_factorial(lint n, lint k);
};

class SpecialBigraph {
public:

    //析构会出问题，所以加一个tag
    //根据是否prepare析构
    bool if_prepare = false;
    
    unsigned num_edges;
    unsigned num_vertices;
    
    unsigned n_vertices[2];
    map<unsigned, unsigned> vertices[2];
    
    vector<unsigned> vertices_in_left;
    vector<unsigned> vertices_in_right;
    
    vector<unsigned> all_vertices;
    
    //unordered_map<int, int> vertexranks;
    unordered_map<int, int> vertexids;
    
    unsigned *deg;
    vector< vector <unsigned> > adj_vec;
    //vector< vector <unsigned> > two_hop_adj_vec;
    unsigned **two_hop_adj_vec;
    unsigned *two_hop_adj_size;
    unsigned *two_hop_adj_maxsize;
    set< pair <unsigned, unsigned> > edges;
    vector< pair <unsigned, unsigned> > list_of_edges;
    
    vector< pair<vector<unsigned>, vector<unsigned> > > bcliques;
    lint largest_index_in_partition[2];
    
    unsigned p;
    unsigned q;
    unsigned n;
    unsigned priority;
    
    unsigned *d;
    unsigned *cd;
    unsigned *adj;
    unsigned **two_hop_d;
    unsigned *two_hop_cd;
    unsigned *two_hop_adj;
    unsigned two_hop_core;
    unsigned core;
    lint result_count;
    
    //static unsigned* old_id;
    //static unsigned* new_id;
    
    class Subgraph{
    public:
        unsigned *n;  // n[l]: number of nodes in two_hop_graph[l]
        unsigned **two_hop_d;
        unsigned *two_hop_adj;
        unsigned *lab;
        unsigned **nodes; //sub[l]: nodes in two_hop_graph[l]
        unsigned **op_vertices;
        unsigned *op_size;
        unsigned two_hop_core;
        
        unsigned *old_id;

        #ifdef FREESUB
        unsigned freel;
        #endif

        
        Subgraph(){}
        ~Subgraph(){

            #ifdef FREESUB
            
            delete[] n;
            for(unsigned i = 1; i <= freel; i++){
                delete[] two_hop_d[i];
                delete[] nodes[i];
                delete[] op_vertices[i];
            }
            delete[] two_hop_d;
            // delete[] two_hop_cd;
            delete[] two_hop_adj;
            delete[] lab;
            delete[] nodes;
            delete[] op_vertices;
            delete[] op_size;

            #endif
            

        }
    };
    
private:
    string data_file_path;
    
    bool anchor_left;
    
    void clear_everything();
    bool all_num(string &s);
    void add_vertex(unsigned A, unsigned side);
    void add_edge(unsigned &A, unsigned &B);
    void get_index(unsigned &A, unsigned side);
    
    void find_combinations(vector< vector<unsigned> > &combs, lint &comb_count, vector<unsigned> &seed_vertices, vector<unsigned> &comb, unsigned beg_offset, unsigned curr_depth, unsigned total_depth);
    void reformat_graph();
    
    void collect_two_hop_adj();
    
    void trim_graph_by_core();
    
    void trim_graph_by_two_hop();
    
    //lint cost_model();
    
    double estimate_cost(unsigned side);
    
    void sort_vertices(unsigned strategy);
    
    void pqclique_thread(unsigned l, Subgraph* sg, lint* rc);
    
    void pqclique_main(unsigned l);
    
    Subgraph* allocsub(unsigned l);
    
    void mksub(unsigned u, Subgraph* sg, unsigned l);
    
    // tools
    void print(vector<unsigned> vec);
    void print_adj();
    void print_adj(unsigned i);
    void print_two_hop_adj();
    void print_edges();
    void print_deg();
    void print_bclique(vector<unsigned> &left_vertices, vector<unsigned> &right_vertices);
    
public:

    SpecialBigraph(string path, unsigned p_value, unsigned q_value, unsigned priority_strategy);

    ~SpecialBigraph();

    void read_graph();

    bool prepare_graph();
    
    void listing_cliques();
    
    void print_results();

//---------------------------------------------------------------------

    void read_graph(vector<pair<uint, uint>> & edgeG);

    SpecialBigraph(unsigned p_value, unsigned q_value, unsigned priority_strategy);
    
};

long long  getCountResult(string path , int p , int q , int strategy, int threadN);

long long  getCountResult(vector<pair<uint, uint>> & edges  , int p , int q , int strategy, int threadN);

#endif
