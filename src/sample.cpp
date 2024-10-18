#include "graph.h"
#include "sample.h"
#include <algorithm>
#include <boost/math/cstdfloat/cstdfloat_types.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/multiprecision/detail/empty_value.hpp>
#include <boost/multiprecision/detail/integer_ops.hpp>
#include <boost/multiprecision/fwd.hpp>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <mutex>
#include <numeric>
#include <ostream>
#include <pthread.h>
#include <queue>
#include <math.h>
#include <random>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp> 
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <shared_mutex>
#include <thread>
#include <unistd.h>
#include <vector>
#include <stack>
#include <numeric>
#define  __SampleTimes 10000000;


#ifdef BOUND_COUNT 
namespace boundCount {
    std::mutex mtxQEqual2;
    size_t QEqual2 = 0;
    std::mutex mtxQEqual1;
    size_t QEqual1 = 0;
    std::mutex mtxHitBound;
    size_t HitBound = 0;
}
#endif

#ifdef SAMPLE_HIT
    std::mutex sampleHit;
    size_t hitTimes = 0;
    std::mutex sampleTime;
    size_t allTime = 0;
#endif




ShadowFinderReturnINFOWithoutE ShadowFinder(tools::CGraph * GlobalG ,tools::Shadow * OriginS,VertexIdx InduceNode){
    

    if(GlobalG->GetDegree(InduceNode)<OriginS->p) return {-1,NULL};
    
    VertexIdx vfirst , vsecond ;
    tools::Shadow * InduceG = tools::ShadowInduceByNodeR(GlobalG, OriginS, InduceNode,vfirst , vsecond);


    if (InduceG == NULL) return {-1,NULL};

    if (InduceG->nVerticesL<InduceG->p or InduceG->nVerticesR<InduceG->q){
        delete InduceG;
        return {-1,NULL};
    }


    if(InduceG->q == 1){

        VertexIdx begin = 0,end = 0;
        auto it1 = InduceG->VerticesR;

        for(; it1 < InduceG->VerticesR+InduceG->nVerticesR ; it1++ ){

            GlobalG->GetNebOfNode(*it1, begin, end);

            if(end - begin < OriginS->q) continue;

            int DegreeOfNodeInInduceG = 0 ;
            for(;begin<end;begin++){
                if(tools::BinarySearch(InduceG->VerticesL, InduceG->nVerticesL , GlobalG->col[begin]) != -1) {
                    DegreeOfNodeInInduceG++;
                    if(DegreeOfNodeInInduceG==InduceG->p) {

                        #ifdef  BOUND_COUNT
                        boundCount::mtxQEqual1.lock();
                        boundCount::QEqual1++;
                        boundCount::mtxQEqual1.unlock();
                        #endif

                        return {1,InduceG};
                    }
                }
            } 
        }

        delete InduceG;
        return {-1,NULL};

    }


    if (InduceG->IfShadow()) 
    { 
        #ifdef BOUND_COUNT   
        boundCount::mtxHitBound.lock();
        boundCount::HitBound++;
        boundCount::mtxHitBound.unlock();
        #endif

        return {1,InduceG};
    }
    else{

        if(InduceG->q == 2){

            VertexIdx begin1,end1,begin2,end2;

            GlobalG->GetNebOfNode(vfirst, begin1, end1);
            GlobalG->GetNebOfNode(vsecond, begin2, end2);

            auto findR = [&GlobalG , &InduceG ,&begin2 ,&end2](VertexIdx id)->bool{
                using namespace tools; 
                return  BinarySearch(GlobalG->col+begin2, end2 - begin2  , id) != -1 
                        && BinarySearch(InduceG->VerticesL,InduceG->nVerticesL , id) != -1
                        ? 1 : 0;
            };

            int NE  =  0 ;

            for(VertexIdx _ = begin1 ; _ < end1 ; _++){
                if(findR(GlobalG->col[_]))  NE++;
                if(NE == InduceG->p){
                    #ifdef BOUND_COUNT
                    boundCount::mtxQEqual2.lock();
                    boundCount::QEqual2++;
                    boundCount::mtxQEqual2.unlock();
                    #endif
                    return {1,InduceG};                    
                }
            }
        }

        return {0,InduceG};

    }        

    return {-2,NULL};

}

bool tools::Shadow::IfShadow(){

    auto getLimit = [](VertexIdx m, VertexIdx n, VertexIdx s , VertexIdx t){
        double limit = std::numeric_limits<double>::max();
        for(int k = 0 ; k <= (s-2) ; k++){
            limit = std::min(limit , pow(s-k-1,1.0/t)*n*pow(m,1.0-1.0/t)  +k*n + (t-1)*pow(m,1+double(k)/t) );
        }
        return limit;
    };

    double limit = getLimit(this->nVerticesL, this->nVerticesR, this->p, this->q);
    return (double)nEdges > limit ; 
    // return false;

}

double tools::Shadow::getLimit(){
    auto getLimit = [](VertexIdx m, VertexIdx n, VertexIdx s , VertexIdx t){
        double limit = std::numeric_limits<double>::max();
        for(int k = 0 ; k <= (s-2) ; k++){
            limit = std::min(limit , pow(s-k-1,1.0/t)*n*pow(m,1.0-1.0/t)  +k*n + (t-1)*pow(m,1+double(k)/t) );
        }
        return limit;
    };
    double limit = getLimit(this->nVerticesL, this->nVerticesR, this->p, this->q); 
    return limit;
}


int intersection(int * vec1, int * vec2, int offset1, int offset2){
    // we assume vec1 and vec2 are sorted by ascending order of the elements
    unsigned i, j;
    i = 0;
    j = 0;
    int re = 0 ;
    while(i < offset1 && j < offset2){
        if(vec1[i] == vec2[j]){
            re ++ ;
            i++;
            j++;
        }
        else if(vec1[i] < vec2[j]){
            i++;
        }
        else{
            j++;
        }
    }
    return re;
}

int countSameElements(const VertexIdx* nums1, int len1, const VertexIdx* nums2, int len2) {

    int count = 0;
    VertexIdx ptr1 = 0, ptr2 = 0;

    while (ptr1 < len1 && ptr2 < len2) {
        if (nums1[ptr1] == nums2[ptr2]) {
            count++;
            ptr1++;
            ptr2++;
        } else if (nums1[ptr1] < nums2[ptr2]) {
            ptr1++;
        } else {
            ptr2++;
        }
    }

    return count;

}

bool IfConstructedClique(tools::CGraph * OriginG , VertexIdx * L , VertexIdx * R , int P , int Q){
    int c = 0 ;
    int c2 = 0;

    int begin , end;
    // std::sort(L,L+P);
    std::sort(R,R+Q);

    for(int _ = 0 ; _ < P ; _++){
        OriginG->GetNebOfNode(*(L+_), begin, end);
        // assert(intersection(R, OriginG->col+begin,Q, end-begin) == countSameElements(R, Q,OriginG->col+begin, end - begin));
        c += countSameElements(R, Q,OriginG->col+begin, end - begin) ;
    }

    #ifdef SAMPLE_HIT
    
    sampleTime.lock();
    allTime++;
    sampleTime.unlock();

    if(c == P*Q){
        sampleHit.lock();
        hitTimes++;
        sampleHit.unlock();
    }
    #endif

    return c == P*Q;
}




