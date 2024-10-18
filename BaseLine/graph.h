

#ifndef GRAPH 
#define GRAPH

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

#define COUNT_ONLY

//#define IS_DEBUGGING

#define TRIM_BY_TWO_HOP_NEIGHBORS

#define USE_CORE_REDUCTION



#define SZ(x) ((int)x.size())
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))


using namespace std;

typedef unsigned long long lint;

class graph{
public:
    unsigned num_edges;
    unsigned num_vertices;
    
    unsigned n_vertices[2];
    map<unsigned, unsigned> vertices[2];
    
    vector<unsigned> vertices_in_left;
    vector<unsigned> vertices_in_right;
    
    vector<unsigned> all_vertices;
    
    unordered_map<int, int> vertexids;
    
    unsigned *deg;
    vector< vector <unsigned> > adj_vec;
    set< pair <unsigned, unsigned> > edges;
    vector< pair <unsigned, unsigned> > list_of_edges;
    
    lint largest_index_in_partition[2];
    
    unsigned p;
    unsigned q;
    unsigned n;
    unsigned priority;
    
    
    
private:

    string data_file_path;
    
    
    void clear_everything();
    bool all_num(string &s);
    void add_vertex(unsigned A, unsigned side);
    void add_edge(unsigned &A, unsigned &B);
    void get_index(unsigned &A, unsigned side);
    

    
    //lint cost_model();
    
    
    // tools
    void print(vector<unsigned> vec);
    void print_adj();
    void print_adj(unsigned i);
    void print_two_hop_adj();
    void print_edges();
    void print_deg();
    void print_bclique(vector<unsigned> &left_vertices, vector<unsigned> &right_vertices);
    
public:

    graph(string path, unsigned p_value, unsigned q_value, unsigned priority_strategy);

    ~graph();

    void read_graph();

    bool prepare_graph();
    
    void listing_cliques();
    
    void print_results();

    void induceSubGraph(vector<pair<uint, uint>>& reEdges,uint nR);
    void induceSubGraph(vector<pair<uint,uint>>& reEdges,uint eL , uint eR);
    void Espar(vector<pair<uint, uint>>& reEdges,double p);
    void CLRSpar(vector<pair<uint, uint>>& reEdges, vector<int>& color);


    
};

#endif
