#include "bclist_parallel.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <random>
#include <sys/types.h>
#include <unordered_map>
#include <utility>
#include <vector>
#include "graph.h"
#include "omp.h"


uint SAMPLETIME = 100000;


struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ (hash2 << 1);
    }
};

double  sampleByEdge(graph& a,int argc,char * argv[]){

    random_device rd;
    mt19937 rand(rd());
    uniform_int_distribution<> distribute(0,a.list_of_edges.size()-1);
    #define genOne (distribute(rand))

    double all = 0;
    all = 0 ;
    long long count = 0;

    unordered_map<pair<uint, uint>, long long, pair_hash> cache;
    
    for(int i = 0 ; i < SAMPLETIME ; i++){
        int nodeInduce = genOne;
        if(cache.find({a.list_of_edges[nodeInduce].first,a.list_of_edges[nodeInduce].second})==cache.end()){
                vector<pair<uint, uint>> subgraph;
            a.induceSubGraph(subgraph,a.list_of_edges[nodeInduce].first,a.list_of_edges[nodeInduce].second);
            count+=cache[{a.list_of_edges[nodeInduce].first,a.list_of_edges[nodeInduce].second}]=\
                getCountResult(subgraph,atoi(argv[2])-1 , atoi(argv[3])-1 ,atoi(argv[4]),atoi(argv[5]));

        }
        all+=cache[{a.list_of_edges[nodeInduce].first,a.list_of_edges[nodeInduce].second}];
        
    }
    cout<<"count："<<count<<endl;
    return (double)a.num_edges*((double)all / (double) (SAMPLETIME));

    
};

double sampleByNode(graph& a,int argc,char * argv[]){

    random_device rd;
    mt19937 rand(rd());
    uniform_int_distribution<> distribute(a.n_vertices[0],a.n_vertices[0]+a.n_vertices[1]-1);
    #define genOne (distribute(rand))
    
    vector<pair<uint, uint>> subgraph;


    double all = 0;
    all = 0 ;

    unordered_map<uint, long long > cahce;

    for(int i = 0 ; i < SAMPLETIME ; i++){

        uint nodeInduce = genOne;
        subgraph.resize(0);
        if(cahce.find(nodeInduce) == cahce.end()){
            a.induceSubGraph(subgraph,nodeInduce);
            cahce[nodeInduce]=getCountResult(subgraph,atoi(argv[2]) , atoi(argv[3])-1 ,atoi(argv[4]),atoi(argv[5]));
        }
        all+=cahce[nodeInduce];

    }

    return (double)a.n_vertices[1]*(all/(double)SAMPLETIME);
    
};

double Espar(graph& a,int argc,char * argv[]){
    vector<pair<uint, uint>> subgraph;
    subgraph.reserve(a.num_edges);
    
    double p = 0.1;
    double all = 0;

    for(int i = 0 ; i < SAMPLETIME ; i++){
        subgraph.resize(0);
        a.Espar(subgraph,p);
        all+=getCountResult(subgraph,atoi(argv[2]) , atoi(argv[3]) ,atoi(argv[4]),atoi(argv[5]));

    }

    return (all/(double)SAMPLETIME)/(pow(p, atoi(argv[2])*atoi(argv[3]))); 

}

double CLRSpar(graph& a,int argc,char * argv[]){

    vector<pair<uint, uint>> subgraph;
    subgraph.reserve(a.num_edges);
        
    int N = 3;
    int SAMPLETIME = 1000000;


    vector<int> color(a.n_vertices[0]+a.n_vertices[1]);


    random_device rd;
    mt19937 rand(rd());
    uniform_int_distribution<> distribute(0,N-1);
    #define genColor (distribute(rand)) 


    double all = 0;



    for(int i = 0 ; i < SAMPLETIME ; i++){
    
        for(auto & n : color){
            n=genColor;
        }

        subgraph.resize(0);
        a.CLRSpar(subgraph, color); ;
        all+=getCountResult(subgraph,atoi(argv[2]) , atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));
   
    }
    
    double p = 1.0/N;

    return ((all)/SAMPLETIME)/pow(p,(atoi(argv[2])*atoi(argv[3]) -1 )); 


}

int main(int argc,char * argv[]){

    graph a(argv[1],atoi(argv[2]) , atoi(argv[3]) ,atoi(argv[4]) );
    // path p q reorder ways
    if(argc<7){
        printf("check input\n");
        exit(1);
    }

    a.read_graph();

    double re ;
    double beg_time,end_time,elapsed_time;

    beg_time = time(NULL);

    switch (atoi(argv[6])) {
        case 2: 
            re= Espar(a, argc, argv);
            break;
        case 1:
            re= sampleByEdge(a, argc, argv);
        case 0:
            re= sampleByNode(a, argc, argv);
    }

    end_time = time(NULL);
    elapsed_time = (end_time - beg_time);
    // cout<<sampleByNode(a, argc, argv);
    printf("time:%lf\n",elapsed_time);
    printf("count:%lf\n",re);    

}



