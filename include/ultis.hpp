

#ifndef bclist_parallel_h
#define bclist_parallel_h

#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <time.h>
#include <cmath>
#include "graph.h"

using namespace std;

typedef unsigned long long lint;

class SpecialBigraph {

    public:
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
    
    unsigned *adj;
    unsigned core;

    vector<unsigned> sortImp;
    
    private:

    string data_file_path;
    
    
    void clear_everything();
    bool all_num(string &s);
    void add_vertex(unsigned A, unsigned side);
    void add_edge(unsigned &A, unsigned &B);
    void get_index(unsigned &A, unsigned side);
    
    void reformat_graph();
    
    void trim_graph_by_core();
    
    void sort_vertices(unsigned strategy);
    void sort_vertices_v2();
    
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
    bool prepare_graph_lv2();
    bool onlysort_graph();
    bool sortGraphWithRecord();
    tools::CGraph * sbgraphToShadow();
    
    
};




#endif 
