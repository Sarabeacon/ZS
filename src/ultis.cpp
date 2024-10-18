
#include "ultis.hpp"
#include "graph.h"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <utility>
#include <vector>


#define SZ(x) ((int)x.size())
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define COUNT_ONLY
#define TRIM_BY_TWO_HOP_NEIGHBORS
#define USE_CORE_REDUCTION


SpecialBigraph::SpecialBigraph(string path, unsigned p_value, unsigned q_value, unsigned priority_strategy){
    clear_everything();
    data_file_path = path;
    p = p_value;
    q = q_value;
    priority = priority_strategy;
}

void SpecialBigraph::clear_everything(){
    largest_index_in_partition[0] = largest_index_in_partition[1] = 0;
    num_vertices = 0;
    num_edges = 0;
    bcliques.clear();
    edges.clear();
    vertices[0].clear(); vertices[1].clear();
    vertices_in_left.clear();
    vertices_in_right.clear();
    adj_vec.clear();
}

bool SpecialBigraph::all_num(string &s) {
    for (int i = 0; i < SZ(s); i++) if ((s[i] >= '0' && s [i] <= '9') == false) return false;
    return true;
}

void SpecialBigraph::add_vertex(unsigned A, unsigned side){
    if(vertices[side].find(A) == vertices[side].end()){
        if(side == 0)
            vertices_in_left.emplace_back(A);
        else
            vertices_in_right.emplace_back(A);
        vertices[side][A] = 1;
    }
}

void SpecialBigraph::add_edge(unsigned &A, unsigned &B){
    add_vertex(A, 0);
    add_vertex(B, 1);
    if(edges.find(make_pair(A, B)) == edges.end()){
        edges.insert(make_pair(A, B));
        num_edges++;
    }
}

void SpecialBigraph::get_index(unsigned &A, unsigned side){
    if(vertices[side].find(A) == vertices[side].end()){
        vertices[side][A] = largest_index_in_partition[side]++;
    }
    A = vertices[side][A];
}

SpecialBigraph::~SpecialBigraph(){
    delete[] deg;
    vertices[0].clear();
    vertices[1].clear();
    vertices_in_left.clear();
    vertices_in_right.clear();
    for(unsigned i = 0; i < adj_vec.size(); i++){
        adj_vec[i].clear();
    }
    adj_vec.clear();
    edges.clear();
    list_of_edges.clear();
    

}

void SpecialBigraph::reformat_graph(){
    // re-format the bipartite graph since we might have removed some vertices
    edges.clear();
    for(unsigned i = 0; i < n_vertices[0]; i++){
        for(unsigned j = 0; j < deg[i]; j++){
            edges.insert(make_pair(i, adj_vec[i][j]));
        }
    }
    num_edges = edges.size();
    
    unsigned removed_count[2];
    removed_count[0] = removed_count[1] = 0;
    unsigned side;
    for(unsigned i = 0; i < num_vertices; i++){
        if(deg[i] == 0){
            side = i < n_vertices[0] ? 0 : 1;
            removed_count[side]++;
        }
    }
    
    n_vertices[0] = n_vertices[0] - removed_count[0];
    n_vertices[1] = n_vertices[1] - removed_count[1];
    
    num_vertices = n_vertices[0] + n_vertices[1];
    
    list_of_edges.clear();
    for(unsigned i = 0; i < adj_vec.size(); i++){
        adj_vec[i].clear();
    }
    adj_vec.clear();
    
    //cout << num_vertices << endl;
    
    vertices[0].clear();
    vertices[1].clear();
    delete[] deg;
    deg = new unsigned[num_vertices];
    memset(deg, 0, num_vertices * sizeof(unsigned));
    adj_vec.resize(num_vertices, vector<unsigned>());
    
    largest_index_in_partition[0] = 0;
    largest_index_in_partition[1] = n_vertices[0];
    for(set < pair <unsigned, unsigned> >::iterator it = edges.begin(); it != edges.end(); it++){
        
        unsigned A = it->first;
        unsigned B = it->second;
        
        //cout << A << " " << B << endl;
        
        get_index(A, 0);
        get_index(B, 1);
        
        //cout << "new : " << A << " " << B << endl;
        
        adj_vec[A].emplace_back(B);
        adj_vec[B].emplace_back(A);
        list_of_edges.emplace_back(make_pair(A, B));
        deg[A]++;
        deg[B]++;
    }
    
    for(unsigned i = 0; i < num_vertices; i++){
        sort(adj_vec[i].begin(), adj_vec[i].end());
    }
    
    #ifdef IS_DEBUGGING
    print_adj();
    
    print_edges();
    #endif
    
    printf("end of reformat graph\n");
    
    printf("#vertices = %u, #left_vertices=%u, #right_vertices=%u\n",num_vertices,n_vertices[0],n_vertices[1]);
    printf("#edges = %lu\n",edges.size());
}

void SpecialBigraph::read_graph(){

    std::ifstream aaa;
    aaa.open(data_file_path.c_str(),std::ios::in);
    assert(aaa.is_open());
    string s;
    while(getline(aaa, s)){
        stringstream ss;
        ss << s;
        vector <string> vec_str;
        for(string z; ss >> z; vec_str.push_back(z));
        if(SZ(vec_str) >= 2){
            bool is_all_number = true;
            for(unsigned i = 0; i < MIN(2, SZ(vec_str)); i++)
                is_all_number &= all_num(vec_str[i]);
            if(is_all_number){
                unsigned A, B;
                ss.clear(); ss << vec_str[0]; ss >> A;
                ss.clear(); ss << vec_str[1]; ss >> B;
                add_edge(A, B);
            }
        }
    }
    vertices[0].clear();
    vertices[1].clear();
    
    n_vertices[0] = vertices_in_left.size();
    n_vertices[1] = vertices_in_right.size();
    num_vertices = n_vertices[0] + n_vertices[1];
    num_edges = edges.size();
    deg = new unsigned[num_vertices]();
    adj_vec.resize(num_vertices, vector<unsigned>());
    
    largest_index_in_partition[0] = 0;
    largest_index_in_partition[1] = n_vertices[0];
    for(set < pair <unsigned, unsigned> >::iterator it = edges.begin(); it != edges.end(); it++){
        unsigned A = it->first;
        unsigned B = it->second;
        get_index(A, 0);
        get_index(B, 1);
        adj_vec[A].emplace_back(B);
        adj_vec[B].emplace_back(A);
        list_of_edges.emplace_back(make_pair(A, B));
        deg[A]++;
        deg[B]++;
    }
    
    #ifdef IS_DEBUGGING
    //print_adj();
    
    //print_edges();
    #endif
    
    printf("end of read graph\n");
    printf("#vertices = %u, #left_vertices = %u, #right_vertices = %u",num_vertices,n_vertices[0],n_vertices[1]);
    printf("#edges = %d\n",num_edges);
}




void SpecialBigraph::trim_graph_by_core(){
    //cout << p << " " << q << endl;
    
    unsigned start_idx, end_idx;
    start_idx = end_idx = 0;
    
    unsigned *to_remove_vertices = new unsigned[num_vertices];
    bool *removed = new bool[num_vertices];
    
    for(unsigned i = 0; i < num_vertices; i++){
        to_remove_vertices[i] = 0;
        removed[i] = false;
    }
    
    // step 1: collect the initial vertices with degree less than q in left and p in right
    for(unsigned i = 0; i < num_vertices; i++){
        if((i < n_vertices[0] && deg[i] < q) || (i >= n_vertices[0] && deg[i] < p)){
            to_remove_vertices[end_idx++] = i;
            removed[i] = true;
        }
    }
        
    // step 2: recursively remove all vertices with degree less than q in left and p in right, i.e., (q,p)-core
    while(start_idx != end_idx){
        unsigned vertex = to_remove_vertices[start_idx++];
        
        //cout << "remove : " << vertex << endl;
        for(int i = 0; i < adj_vec[vertex].size(); i++){
            unsigned other = adj_vec[vertex][i];
            if(!removed[other]){
                deg[other]--;
                
                if((other < n_vertices[0] && deg[other] < q) || (other >= n_vertices[0] && deg[other] < p)){
                    to_remove_vertices[end_idx++] = other;
                    removed[other] = true;
                }
            }
        }
        adj_vec[vertex].clear();
        deg[vertex] = 0;
    }
        
    for(unsigned i = 0; i < num_vertices; i++){
        if(!removed[i]){
            vector<unsigned> new_vec(deg[i]);
            unsigned idx = 0;
            for(unsigned j = 0; j < adj_vec[i].size(); j++){
                if(!removed[adj_vec[i][j]]){
                    new_vec[idx++] = adj_vec[i][j];
                }
            }
            adj_vec[i].clear();
            adj_vec[i] = new_vec;
        }
    }
    
    #ifdef IS_DEBUGGING
    print_adj();
    #endif
    
    printf("after trimming by core\n");
    
    reformat_graph();

    printf("end of trim graph by core\n");
}

void SpecialBigraph::sort_vertices(unsigned strategy = 0){

    switch (strategy) {

        case 0: {

            vector<pair<unsigned,unsigned>>  order;
            for(unsigned A = n_vertices[0] ; A < num_vertices ; A++)
                order.push_back(make_pair(A, deg[A]));
            


            auto sortByDegree = [](const pair<unsigned, unsigned>& A , pair<unsigned, unsigned>& B){return A.second < B.second;};
            
            std::sort(order.begin(),order.end(),sortByDegree);

            vector<unsigned> impR(n_vertices[1]);

            for(unsigned i = 0 ; i < n_vertices[1] ; i ++){
                impR[order[i].first - n_vertices[0]] = i;
            }
            
                unsigned oldId = 100;
                unsigned newId = impR[oldId];

            
            edges.clear();

            for(unsigned i = 0 ; i < n_vertices[0] ; i ++){
                for(unsigned j : adj_vec[i] ){
                    edges.insert(make_pair(i, impR[j-n_vertices[0]] + n_vertices[0]));
                }
            }

            num_edges = num_edges+0;
            n_vertices[0] = n_vertices[0]+0;
            n_vertices[1] = n_vertices[1]+0;
            num_vertices = n_vertices[0] + n_vertices[1];
            
            list_of_edges.clear();
            for(unsigned i = 0; i < adj_vec.size(); i++){
                adj_vec[i].clear();
            }
            adj_vec.clear();
            adj_vec.resize(num_vertices, vector<unsigned>());
            
            for(set < pair <unsigned, unsigned> >::iterator it = edges.begin(); it != edges.end(); it++){
                unsigned A = it->first;
                unsigned B = it->second;
                adj_vec[A].emplace_back(B);
                adj_vec[B].emplace_back(A);
                list_of_edges.emplace_back(make_pair(A, B));
                deg[A]++;
                deg[B]++;
            }
            
            for(unsigned i = 0; i < num_vertices; i++){
                sort(adj_vec[i].begin(), adj_vec[i].end());
            }
            
            
            printf("end of resort graph , resort by degree in nodeR\n");
            printf("#vertices = %u, #left_vertices=%u, #right_vertices=%u\n",num_vertices,n_vertices[0],n_vertices[1]);
            printf("#edges = %lu\n",edges.size());

        }
        case 1:{
            
        }
        default:
            break;
    }

    

}

bool SpecialBigraph::prepare_graph(){

    trim_graph_by_core();
    
    if(n_vertices[0] == 0 || n_vertices[1] == 0){
        printf("No results because the graph is pruned by core\n");
        return false;
    }

    sort_vertices(priority);
    // sort_vertices(1);
    
    printf("finish sorting vertices\n");
    
    return true;
}

bool SpecialBigraph::onlysort_graph(){

    sort_vertices(priority);

    printf("only sort graph\n");
    
    return true;
}

tools::CGraph *   SpecialBigraph::sbgraphToShadow(){

    tools::CGraph * re = new tools::CGraph;

    re->p = this->p;
    re->q = this->q;
    re->nEdges = this->num_edges;
    re->nVerticesL = this->n_vertices[0];
    re->nVerticesR = this->n_vertices[1];

    re->ptr = new VertexIdx[this->num_vertices+1]; 
    re->col = new VertexIdx[this->num_edges * 2];    

    re->ptr[0] = 0 ;

    unsigned ptrIndex = 1;
    unsigned colIndex = 0;
    
    for(auto i: this->adj_vec){
        re->ptr[ptrIndex] = re->ptr[ptrIndex-1] + i.size();
        ptrIndex++;
        for(auto b : i){
            re->col[colIndex] = b;
            colIndex++;
        }
    }

    return re;

}


bool SpecialBigraph::prepare_graph_lv2(){

    trim_graph_by_core();
    
    if(n_vertices[0] == 0 || n_vertices[1] == 0){
        printf("No results because the graph is pruned by core\n");
        return false;
    }
    
    //abcore


    sort_vertices(priority);
    
    printf("finish sorting vertices\n");
    
    return true;
}

void SpecialBigraph::sort_vertices_v2(){

    int strategy = 0;
    switch (strategy) {

        case 0: {

            vector<pair<unsigned,unsigned>>  order;
            for(unsigned A = n_vertices[0] ; A < num_vertices ; A++)
                order.push_back(make_pair(A, deg[A]));
            
            // for(auto A : order){
            //     assert(A.first >= n_vertices[0]);
            //     // assert(!(A.first == 709));
            // }

            auto sortByDegree = [](const pair<unsigned, unsigned>& A , pair<unsigned, unsigned>& B){return A.second < B.second;};
            
            std::sort(order.begin(),order.end(),sortByDegree);

            vector<unsigned> impR(n_vertices[1]);

            for(unsigned i = 0 ; i < n_vertices[1] ; i ++){
                impR[order[i].first - n_vertices[0]] = i;
                this->sortImp[order[i].first] = i+ this->n_vertices[0];
            }
            
                unsigned oldId = 100;
                unsigned newId = impR[oldId];
            //

            
            edges.clear();

            for(unsigned i = 0 ; i < n_vertices[0] ; i ++){
                for(unsigned j : adj_vec[i] ){
                    edges.insert(make_pair(i, impR[j-n_vertices[0]] + n_vertices[0]));
                }
            }

            num_edges = num_edges+0;
            n_vertices[0] = n_vertices[0]+0;
            n_vertices[1] = n_vertices[1]+0;
            num_vertices = n_vertices[0] + n_vertices[1];
            
            list_of_edges.clear();
            for(unsigned i = 0; i < adj_vec.size(); i++){
                adj_vec[i].clear();
            }
            adj_vec.clear();
            adj_vec.resize(num_vertices, vector<unsigned>());
            
            for(set < pair <unsigned, unsigned> >::iterator it = edges.begin(); it != edges.end(); it++){
                unsigned A = it->first;
                unsigned B = it->second;
                adj_vec[A].emplace_back(B);
                adj_vec[B].emplace_back(A);
                list_of_edges.emplace_back(make_pair(A, B));
                deg[A]++;
                deg[B]++;
            }
            
            for(unsigned i = 0; i < num_vertices; i++){
                sort(adj_vec[i].begin(), adj_vec[i].end());
            }
            
            
            printf("end of resort graph , resort by degree in nodeR\n");
            printf("#vertices = %u, #left_vertices=%u, #right_vertices=%u\n",num_vertices,n_vertices[0],n_vertices[1]);
            printf("#edges = %lu\n",edges.size());

        }
        default:
            break;
    }

    


}

bool SpecialBigraph::sortGraphWithRecord(){


    this->sortImp.resize(this->num_vertices);
    this->sort_vertices_v2();


}

