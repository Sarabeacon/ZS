#ifndef SHADOW 
#define SHADOW

#include "./graph.h"
#include <algorithm>
#include <boost/multiprecision/cpp_int/cpp_int_config.hpp>
#include <cassert>
#include <cstddef>
#include <list>
#include <mutex>
#include <numeric>
#include <random>
#include <vector>
#include "setting.h"



//计数工具
#ifdef BOUND_COUNT 
namespace boundCount {
    extern std::mutex mtxQEqual2;
    extern size_t QEqual2;
    extern std::mutex mtxQEqual1;
    extern size_t QEqual1;
    extern std::mutex mtxHitBound;
    extern size_t HitBound;
}
#endif

#ifdef MEMORY_USE
namespace memCount {
    extern size_t memoryMax;
    extern size_t memoryTotal;
}
#endif


struct sampleEngine{

    sampleEngine(std::vector<double> ptr):rand(dRand()),uRand(0.0,1.0){
        nodeDistrubute.resize(ptr.size());
        double sum = std::accumulate(ptr.begin(),ptr.end(),(double)0);
        std::for_each(ptr.begin(),ptr.end(),[&sum](double & a){a=a/sum;});
        std::partial_sum(ptr.begin(),ptr.end(),nodeDistrubute.begin());
    };

    sampleEngine(const sampleEngine&) = delete;
    sampleEngine() = delete;
    sampleEngine(sampleEngine&& __t) =delete;

    VertexIdx operator()(){
        auto it = std::upper_bound(nodeDistrubute.begin(),nodeDistrubute.end(),uRand(rand)); 
        return it - nodeDistrubute.begin();
    }

    private:
    std::vector<double> nodeDistrubute;
    std::random_device dRand;
    std::mt19937 rand;
    std::uniform_real_distribution<double> uRand;

};




struct ShadowFinderReturnINFOWithoutE{
    int tag;//
    tools::Shadow * G;
    // ShadowFinderReturnINFOWithoutE& operator=(ShadowFinderReturnINFO && a);
};

void CoreFunction( tools::CGraph *  OriginG,std::list<tools::Shadow*> & Shadow);

ShadowFinderReturnINFOWithoutE ShadowFinder(tools::CGraph * GlobalG ,tools::Shadow * OriginS,VertexIdx InduceNode);

//
// void ShadowFinder_DFS(tools::CGraph *  GlobalG, std::list<tools::Shadow *> & ShadowSave);

// void ShadowFinder_DFS_Multithreading(tools::CGraph *  GlobalG, std::list<tools::Shadow *> & ShadowSave , int ThreadNum);


// long SampleShadow(tools::CGraph * OriginG , std::list<tools::Shadow *> & SL , int P , int Q  , double Epsilon , double Delta );

// long SampleShadow_Multithreading(tools::CGraph * OriginG , std::list<tools::Shadow *> & SL , int P , int Q  , double Epsilon , double Delta , const int ThreadNum);


bool IfConstructedClique(tools::CGraph * OriginG , VertexIdx * L , VertexIdx * R , int P , int Q);

// void SampleIndex(const std::vector<double>& distribution ,std::vector<VertexIdx> & Missions, unsigned long long Times);

// int SampleNodeToConstructClique( tools :: CGraph * OriginG , tools::Shadow * S);



#endif