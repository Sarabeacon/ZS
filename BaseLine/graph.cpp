


#include "graph.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <omp.h>
#include <random>
#include <sys/types.h>
#include <utility>
#include <vector>

// below are special bigraphs

graph::graph(string path, unsigned p_value, unsigned q_value, unsigned priority_strategy){
    clear_everything();
    data_file_path = path;
    p = p_value;
    q = q_value;
    priority = priority_strategy;
}

void graph::clear_everything(){
    largest_index_in_partition[0] = largest_index_in_partition[1] = 0;
    num_vertices = 0;
    num_edges = 0;
    edges.clear();
    vertices[0].clear(); vertices[1].clear();
    vertices_in_left.clear();
    vertices_in_right.clear();
    adj_vec.clear();
}

bool graph::all_num(string &s) {
    for (int i = 0; i < SZ(s); i++) if ((s[i] >= '0' && s [i] <= '9') == false) return false;
    return true;
}

void graph::add_vertex(unsigned A, unsigned side){
    if(vertices[side].find(A) == vertices[side].end()){
        if(side == 0)
            vertices_in_left.emplace_back(A);
        else
            vertices_in_right.emplace_back(A);
        vertices[side][A] = 1;
    }
}

void graph::add_edge(unsigned &A, unsigned &B){
    add_vertex(A, 0);
    add_vertex(B, 1);
    if(edges.find(make_pair(A, B)) == edges.end()){
        edges.insert(make_pair(A, B));
        num_edges++;
    }
}

void graph::get_index(unsigned &A, unsigned side){
    if(vertices[side].find(A) == vertices[side].end()){
        vertices[side][A] = largest_index_in_partition[side]++;
    }
    A = vertices[side][A];
}

void graph::read_graph(){

    std::ifstream aaa;
    aaa.open(data_file_path.c_str(),std::ios::in);

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
    
    // printf("end of read graph\n");
    // printf("#vertices = %u, #left_vertices = %u, #right_vertices = %u",num_vertices,n_vertices[0],n_vertices[1]);
    // printf("#edges = %d\n",num_edges);
}


graph::~graph(){
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


// the following are tool functions

void graph::print(vector<unsigned> vec){
    for(vector<unsigned>::iterator it = vec.begin(); it != vec.end(); it++){
        printf("%u ",*it);
    }
    printf("\n");
}

void graph::print_adj(){
    for(unsigned i = 0; i < num_vertices; i++){
        cout << "adj[" << i << "] : ";
        printf("adj[%u]:",i);
        for(unsigned j = 0; j < deg[i]; j++){
            printf("%d ",adj_vec[i][j]);
        }
        printf("\n");
    }
}

void graph::print_adj(unsigned i){
    for(unsigned j = 0; j < deg[i]; j++){
        printf("%d ",adj_vec[i][j]);
    }
    printf("\n");
}


void graph::print_edges(){
    cout << "#vertices : " << num_vertices <<endl;
    cout << "#edges : " << num_edges << endl;
    for(vector< pair<unsigned, unsigned> >::iterator it = list_of_edges.begin(); it != list_of_edges.end(); it++){
        cout << it->first << " " << it->second << endl;
    }
}

void graph::print_deg(){
    sort(deg, deg+num_vertices);
    for(unsigned i = num_vertices - 1; i >= num_vertices - 1000; i--){
        cout << deg[i] << endl;
    }
}

void graph::print_bclique(vector<unsigned> &left_vertices, vector<unsigned> &right_vertices){
    for(unsigned i = 0; i < left_vertices.size(); i++){
        printf("%d ",left_vertices[i]);
    }
    cout << "| ";
    for(unsigned i = 0; i < right_vertices.size(); i++){
        printf("%d ",right_vertices[i]);
    }
    printf("\n");
}



#define INDUCE

void graph::induceSubGraph(vector<pair<uint, uint>>& reEdges,uint nR){

    #ifdef INDUCE
    #endif

    for(auto & subU : this->adj_vec[nR]){
    
        for(auto & subV : this->adj_vec[subU]){

            if(subV>nR){
                reEdges.push_back({subU,subV});
            }

        }
    }



}

void graph::induceSubGraph(vector<pair<uint, uint>>& reEdges,uint eL , uint eR){

    vector<uint> candL;

    for(auto& subL :  adj_vec[eR]){
        if(subL>eL){
            candL.push_back(subL);
        }
    }

    for(auto & subR : adj_vec[eL]){

        if(subR>eR){
            for(auto & subL : adj_vec[subR]){
                if(binary_search(candL.begin(),candL.end(),subL)){
                    reEdges.push_back({subL,subR});
                }
            } 
        }

        
    }

}

void graph::Espar(vector<pair<uint, uint>>& reEdges , double p ){

    #ifdef INDUCE
    #endif
    
    random_device rd;
    mt19937 rand(rd());
    uniform_real_distribution<double> dist(0,1);
    #define ifChoose (dist(rand)<p)

    for(auto & e : this->list_of_edges){
        if(ifChoose){
            reEdges.push_back(e);
        }
    }

}
void graph::CLRSpar(vector<pair<uint, uint>>& reEdges, vector<int>& color){
    for(auto & e : this->list_of_edges){
        if(color[e.first]==color[e.second]){
            reEdges.push_back(e);
        }
    }
}