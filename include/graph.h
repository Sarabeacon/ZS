#ifndef GRAPH
#define GRAPH

#include <algorithm>
#include <boost/lexical_cast/try_lexical_convert.hpp>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <assert.h>
#include <ios>
#include <iostream>

#include <algorithm>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <list>
#include <vector>
#include <math.h>

typedef int VertexIdx ;
typedef int EdgeIdx ;
typedef std::unordered_map<int, std::list<int>> NODE_INFO;

namespace tools {




VertexIdx BinarySearch(VertexIdx * array , VertexIdx size , VertexIdx id);

VertexIdx FindFirstGreaterThen(VertexIdx * array , VertexIdx end ,VertexIdx id);

bool NodeOrder(VertexIdx if_first_more_than, VertexIdx second);


bool cmpArray(VertexIdx * A, VertexIdx* B ,VertexIdx size );

    struct Graph
    {
        public:
            VertexIdx * row;
            VertexIdx * col;
            VertexIdx nVerticesL;
            VertexIdx nVerticesR;
            VertexIdx nEdges;
            ~Graph();
            void print();
    };
    
    Graph * ReadGraph(std::string PATH);



    struct CGraph
    {
        public:
        VertexIdx   nVerticesL;
        VertexIdx   nVerticesR;
        EdgeIdx     nEdges;
        VertexIdx  *ptr;   
        VertexIdx  *col;   
        int p;
        int q;

        public:
        inline void       GetNebOfNode( VertexIdx  & node, VertexIdx & start ,VertexIdx & end);
        inline VertexIdx  GetDegree( VertexIdx  &  node);
        inline VertexIdx  get2HopAdjNum(VertexIdx  node);
        inline VertexIdx  getAdjNum(VertexIdx  node);
        ~CGraph();
        void print();
    };

    tools::CGraph * CooToCsr(tools::Graph * OriginGraph) ;

//------------New Born 

    struct Shadow{
        public:
        VertexIdx   nVerticesL;
        VertexIdx   nVerticesR;
        EdgeIdx     nEdges;
        int p;
        int q;
        VertexIdx  *VerticesL;   
        VertexIdx  *VerticesR;  
        ~Shadow(){delete [] VerticesL ;delete []  VerticesR ;};
        public:
        bool IfShadow();
        double getLimit();
        void print();
        bool operator ==(Shadow & b){
            bool tag1 = cmpArray(this->VerticesL, b.VerticesL, this->nVerticesL);
            bool tag2 = cmpArray(this->VerticesR, b.VerticesR, this->nVerticesR);
            return this->nVerticesL==b.nVerticesL && this->nVerticesR == b.nVerticesR && tag1 && tag2;
        }

        

    };

    Shadow* ShadowInduceByNodeR(CGraph *OriginG,Shadow * OriginS, VertexIdx Nodex,VertexIdx & vfirst , VertexIdx & vsecond);
    Shadow* csrToShadow(CGraph * a, int p ,int q);



}

inline VertexIdx tools::CGraph::GetDegree(VertexIdx &node){
    return this->ptr[node+1]-ptr[node]; 
}

inline void tools::CGraph::GetNebOfNode( VertexIdx  & node, VertexIdx &start, VertexIdx &end){
    start = this->ptr[node]; 
    end = this->ptr[node+1]; 
}


inline VertexIdx  tools::CGraph::get2HopAdjNum(VertexIdx  node){
    VertexIdx start,end;
    this->GetNebOfNode(node, start,end);
    VertexIdx num = 0;
    std::vector<bool> tempL(this->nVerticesL+this->nVerticesR,0);
    for(;start < end ; start++){
        VertexIdx nodeU = this->col[start];
        VertexIdx s,e; 
        this->GetNebOfNode(nodeU, s, e);
        for(int i = e-1 ; i>=s ;i--){
            if(*(this->col+i) <= node) break;
            tempL[*(this->col+i)] = 1;
        }

    }
    for (int i=0; i<(this->nVerticesL+this->nVerticesR); i++) {
        if(tempL[i]){
            num++;
        }
    }

    return num;

}

inline VertexIdx tools::CGraph::getAdjNum(VertexIdx  node){
    return this->GetDegree(node);
}



#endif