#include <algorithm>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include "graph.h"
#define __MISSIONSHADOW


bool tools::cmpArray(VertexIdx * A, VertexIdx* B ,VertexIdx size ){
    for(int _  = 0 ; _ <size ;_++){
        if(*(A+_) != *(B+_)){
            return false;
        }
    }
    return true;
};

bool tools::NodeOrder(VertexIdx if_first_more_than, VertexIdx second){
    return if_first_more_than>second;
} 

VertexIdx tools::BinarySearch(VertexIdx *array, VertexIdx size, VertexIdx val){


    VertexIdx low = 0;
    VertexIdx high = size-1;
    VertexIdx mid;

    while (low <= high)
    {
        mid = (low+high)/2;
        if (array[mid] == val)
            return mid;
        if (array[mid] > val)
            high = mid-1;
        if (array[mid] < val)
            low = mid+1;
    }
    return -1;
    
}

VertexIdx tools::FindFirstGreaterThen(VertexIdx *array, VertexIdx end, VertexIdx id){
    VertexIdx low = 0;
    VertexIdx high = end -1;
    VertexIdx result = -1;
    VertexIdx media = 0 ;
    while (low <= high) {
        media = low + (high-low)/2;
        if(array[media] > id){
            result = media;
            high = media -1 ;
        }else {
            low = media + 1;
        }
    }

    return result;
}


void tools::Graph::print(){
    for(int _ = 0 ; _ < this->nEdges ; _ ++){
        std::cout << this->row[_] << "  "<<col[_] << "\n";
    }
}

tools::Graph * tools::ReadGraph(std::string PATH ){

    std::fstream FILE;
    FILE.open(PATH);
    assert(FILE.is_open());

    VertexIdx MaxL=-1,MaxR=-1;
    VertexIdx a , b;
    VertexIdx NumE=0;
    while (FILE>>a>>b) {
        MaxL = std::max(MaxL,a);
        MaxR = std::max(MaxR,b);
        NumE++;
    }


    FILE.close();
    FILE.open(PATH);

    std::vector<std::pair<VertexIdx, VertexIdx >> EdgesLToR;
    EdgesLToR.reserve(NumE);
    std::vector<std::pair<VertexIdx, VertexIdx >> EdgesRToL;
    EdgesRToL.reserve(NumE);

    while (FILE>>a>>b) {
        EdgesLToR.push_back(std::pair<VertexIdx, VertexIdx>(a,b+MaxL));
        EdgesRToL.push_back(std::pair<VertexIdx, VertexIdx>(b+MaxL,a));
    } 

    auto compare = [](const std::pair<VertexIdx, VertexIdx>& a , const std::pair<VertexIdx, VertexIdx>& b ){
        if (a.first==b.first) {
            return a.second<b.second;
        }
        return a.first<b.first;
    };
    std::sort(EdgesLToR.begin(),EdgesLToR.end(),compare);
    std::sort(EdgesRToL.begin(),EdgesRToL.end(),compare);

    std::vector<std::pair<VertexIdx, VertexIdx >> Edges = EdgesLToR;
    Edges.resize(NumE*2);
    std::copy(EdgesRToL.begin(),EdgesRToL.end(),Edges.begin()+NumE);


    
    VertexIdx NumOfL = 1;  
    for(VertexIdx _ = 1 ; _ < NumE ; _++){
        if(EdgesLToR[_].first!=EdgesLToR[_-1].first) NumOfL++;
    }
    VertexIdx NumOfR = 1;
    for(VertexIdx _ = 1 ; _ < NumE ; _++){
        if(EdgesRToL[_].first!=EdgesRToL[_-1].first) NumOfR++;
    }

    VertexIdx * LIMP = new VertexIdx[NumOfL];
    VertexIdx IterIdx = 0;
    LIMP[IterIdx] = EdgesLToR[IterIdx].first;
    for(VertexIdx _ = 1 ; _ < NumE ; _++){
        if(EdgesLToR[_-1].first != EdgesLToR[_].first){
            IterIdx++;
            LIMP[IterIdx] = EdgesLToR[_].first;
        };
    }
    VertexIdx * RIMP = new VertexIdx[NumOfR];
    IterIdx = 0;
    RIMP[IterIdx] = EdgesRToL.begin()->first;
    for(VertexIdx _ = 1 ; _ < NumE ;_++){
        if(EdgesRToL[_-1].first != EdgesRToL[_].first){
            IterIdx++;
            RIMP[IterIdx] = EdgesRToL[_].first;
        }
    }

    auto RNodeGetIMP = [&RIMP,&NumOfL,&NumOfR](VertexIdx nodex)->VertexIdx{return BinarySearch(RIMP,NumOfR,nodex)==-1?-1:NumOfL+BinarySearch(RIMP,NumOfR,nodex);};
    auto LNodeGetIMP = [&LIMP,&NumOfL,&NumOfR](VertexIdx nodex){return BinarySearch(LIMP, NumOfL, nodex);};

    for(int _ = 0 ; _ <NumE ; _ ++){
        Edges[_].first = LNodeGetIMP(Edges[_].first);
        Edges[_].second = RNodeGetIMP(Edges[_].second);
    }
    for(int _ = NumE ; _ <2*NumE ; _ ++){
        Edges[_].first = RNodeGetIMP(Edges[_].first);
        Edges[_].second = LNodeGetIMP(Edges[_].second);
    }



    Graph * ReGraphPtr = new Graph({new VertexIdx[2*NumE],new VertexIdx[2*NumE],NumOfL,NumOfR,NumE});

    auto ColPtr = ReGraphPtr->col;
    auto RowPtr = ReGraphPtr->row;

    for(auto & _ : Edges){
        *ColPtr = _.second;
        *RowPtr = _.first;
        ColPtr++;
        RowPtr++;
    }
    

    delete [] RIMP;
    delete [] LIMP;
    return ReGraphPtr;
}


tools::CGraph * tools::CooToCsr(tools::Graph * OriginGraph){

    tools::CGraph * ReNewCGraph =\
    new tools::CGraph({OriginGraph->nVerticesL,OriginGraph->nVerticesR,OriginGraph->nEdges,\
    new VertexIdx[OriginGraph->nVerticesL+OriginGraph->nVerticesR+1],new VertexIdx[OriginGraph->nEdges*2]});


    std::copy(OriginGraph->col,OriginGraph->col+(OriginGraph->nEdges)*2,ReNewCGraph->col);

    auto itr = ReNewCGraph->ptr;
    *itr=0;
    itr++;
    *itr=1;
    for(int _ = 1;_ < ReNewCGraph->nEdges*2 ;_++){
        if(OriginGraph->row[_]==OriginGraph->row[_-1]){
            (*itr)++;
        }
        else {
            *itr=*(itr-1)+*itr;
            itr++;
            *itr=1;
        }
    }
    *itr=*(itr-1)+*itr;

    return ReNewCGraph;

}

tools::Graph::~Graph(){
    delete [] col;
    delete [] row;
}

void tools::CGraph::print(){
    for(VertexIdx _ = 0 ; _ < this->nVerticesL+this->nVerticesR ; _++){
        std::cout<<"id:"<<_<<" D:"<<this->GetDegree(_)<<" Neb: |";
        VertexIdx begin , end ;
        this->GetNebOfNode(_, begin, end);
        for(int __ = begin ; __ < end ; __++){
            std::cout << this->col[__] <<"|";
        }
        std::cout<<"\n";
    }
}

tools::CGraph::~CGraph(){
    delete [] ptr;
    delete [] col;
}


void tools::Shadow::print(){

    std::cout<<"\n\n---------------------\n";
    printf("L:%d R:%d E:%d p:%d  q:%d \n",this->nVerticesL,this->nVerticesR,this->nEdges,this->p,this->q);
    
    std::cout<<"ifShadow:"<<this->IfShadow()<<'\n';
    std::cout<<"limit:"<< 0.5 * pow(p-q+1,1.0/q) * pow((double)nVerticesR,2-1.0/q) + 0.5 * (q-1) * pow((double)nVerticesR , 2-2.0/q) + 0.5 * (q-2) * nVerticesR;
    std::cout<<"  E:"<<this->nEdges<<'\n';
    std::cout<<"L:|";
    for(int _ = 0 ; _ < this->nVerticesL ; _++){
        std::cout << this->VerticesL[_]<<"|";
    }
    std::cout<<'\n';

    std::cout<<"R:|";
    for(int _ = 0 ; _ < this->nVerticesR ; _++){
        std::cout << this->VerticesR[_]<<"|";
    }
    std::cout<<'\n';

}

tools::Shadow*  tools::csrToShadow(tools::CGraph * OriginG , int p ,int q){
    tools::Shadow * StartS = new tools::Shadow\
    ({OriginG->nVerticesL,OriginG->nVerticesR,OriginG->nEdges,OriginG->p,OriginG->q,\
    new VertexIdx[OriginG->nVerticesL],new VertexIdx[OriginG->nVerticesR]});
    for(VertexIdx _ = 0 ; _ < StartS->nVerticesL ; _ ++){
        StartS->VerticesL[_] = _;    
    }
    for(VertexIdx _ = 0 ; _ < StartS->nVerticesR ; _ ++){
        StartS->VerticesR[_] = StartS->nVerticesL + _;    
    }
    StartS->p = p;
    StartS->q = q;
    return StartS;
}

tools::Shadow * tools::ShadowInduceByNodeR(CGraph *OriginG,Shadow * OriginS, VertexIdx Nodex , VertexIdx & vfirst , VertexIdx & vsecond){//

    if(Nodex < OriginG->nVerticesL) throw("check input");

    VertexIdx begin = 0 ;
    VertexIdx end = 0; 
    VertexIdx NumL =0 ;
    VertexIdx NumE = 0;
    
    VertexIdx*  TempNewShadowVerticesL = new VertexIdx[OriginG->GetDegree(Nodex)];
    OriginG->GetNebOfNode(Nodex,begin, end); 
    auto it3 = TempNewShadowVerticesL ;

    for(;begin<end;begin++){
        if(tools::BinarySearch(OriginS->VerticesL,OriginS->nVerticesL,OriginG->col[begin])!=-1){
            *it3 = OriginG->col[begin];
            it3++;
        }
    } 

    NumL = it3 - TempNewShadowVerticesL;

    #ifdef  __MISSIONSHADOW
    if(NumL < OriginS->p){
        delete [] TempNewShadowVerticesL;
        return NULL;
    }
    #endif

    VertexIdx * NewShadowVerticesL = new VertexIdx[NumL];
    std::copy(TempNewShadowVerticesL,it3,NewShadowVerticesL);


    size_t TempNumR = 0 ; 

    for(VertexIdx _  = 0 ; _ < NumL ; _ ++){
        OriginG->GetNebOfNode(NewShadowVerticesL[_], begin, end);
        TempNumR += end - begin;
    }

    VertexIdx * TempVerticesR = new VertexIdx[TempNumR];
    VertexIdx * it1 = TempVerticesR;

    for(VertexIdx _ = 0 ; _ < NumL ; _ ++){
        OriginG->GetNebOfNode(NewShadowVerticesL[_], begin, end);
        auto new_start  = tools::FindFirstGreaterThen(OriginG->col+begin, end - begin, Nodex);//找到第一个偏移
        if(new_start == -1) continue;
        begin += new_start;

        NumE += end - begin ;
        for(VertexIdx __ = begin ; __ < end ; __ ++){
            *it1 = OriginG->col[__] ;
            it1++;
        }
    }

    #ifdef  __MISSIONSHADOW
    if (TempVerticesR == it1) {
        delete []  TempVerticesR; 
        delete []  NewShadowVerticesL;
        delete [] TempNewShadowVerticesL;
        return NULL;
    }
    #endif


    std::sort(TempVerticesR, it1);

        if(OriginS->q==3){
            std::map<VertexIdx, VertexIdx> m ;
            for(VertexIdx index = 0 ; index < NumE ; index++){
                if(m.find(TempVerticesR[index])==m.end()){
                    m[TempVerticesR[index]] = 1;
                }
                else {
                    m[TempVerticesR[index]] ++;
                }
            }
            VertexIdx dfirst = 0  , dsecond = 0 ;
            VertexIdx nfirst = 0  , nsecond = 0 ;
            for(auto & a : m){
                if (a.second > dfirst) {
                    dsecond = dfirst; nsecond = nfirst;
                    nfirst = a.first ; dfirst = a.second;
                }
                else if (a.second > dsecond) {
                    dsecond = a.second ; nsecond = a.first;
                }
            }
            vfirst = nfirst;
            vsecond = nsecond;
        }

    VertexIdx * it2 = TempVerticesR + 1 ; 
    VertexIdx * f = TempVerticesR + 1;
    for(; f < it1 ; f++){
        if (*f == *(f-1)) {

        } 
        else {       
            *it2 = *f;
            it2++; 
        }
    }



    


    VertexIdx NumR = it2 - TempVerticesR;
    #ifndef __MISSIONSHADOW
        if(NumE==0) NumR =0;
    #endif
    //-----------------

    #ifdef  __MISSIONSHADOW
        if(NumR < OriginS->q -1){
            delete [] TempNewShadowVerticesL;
            delete [] TempVerticesR;
            return NULL;
        }
    #endif

    VertexIdx * NewShadowVerticesR = new VertexIdx[NumR]; 

    std::copy(TempVerticesR , it2 ,NewShadowVerticesR);

    auto aaaa =  new tools::Shadow({NumL,NumR,NumE,OriginS->p,OriginS->q-1,NewShadowVerticesL,NewShadowVerticesR});

    delete [] TempVerticesR;
    delete [] TempNewShadowVerticesL;

    return aaaa;

}

