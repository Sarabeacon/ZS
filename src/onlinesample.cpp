#include "graph.h"
#include "sample.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <ios>
#include <iostream>
#include <iterator>
#include <list>
#include <numeric>
#include <ostream>
#include <random>
#include <stack>
#include <thread>
#include <unistd.h>
#include <utility>
#include <vector>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp> 
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "sample.h"
#include "setting.h"

#ifdef MEMORY_USE
namespace memCount {
    size_t memoryMax = 0;
    size_t memoryTotal = 0;
}
#endif


namespace OnlineMT_Induce {

    void ShadowFinder_DFS_Multithreading_Core(tools::CGraph * GlobalG,tools::Shadow * StartS, std::list<tools::Shadow *> & ShadowSave, const int ThreadId,const int ThreadNum){

        for(VertexIdx ptrNodeR = 0 + ThreadId ; ptrNodeR < StartS->nVerticesR ; ptrNodeR += ThreadNum){

            VertexIdx NodeR = *(StartS->VerticesR + ptrNodeR);

            std::stack<tools::Shadow * > T ; 
        
            auto BS = ShadowFinder(GlobalG,StartS,NodeR);
            if(BS.tag == 1){
                ShadowSave.push_back(BS.G);
                continue;
            }else if (BS.tag == 0) {
                T.push(BS.G);
            }
            else {
                continue;
            }


            while (!T.empty()) {
                
                auto NowG = T.top(); 
                T.pop();

                for(VertexIdx _ = 0 ; _ < NowG->nVerticesR ; _ ++){

                    auto INFO = ShadowFinder(GlobalG, NowG, NowG->VerticesR[_]);

                    if (INFO.tag == 1) {
                        ShadowSave.push_back(INFO.G);
                    }else if (INFO.tag == 0) {
                        T.push(INFO.G);
                    }
                    else {
                        continue;
                    }
                }

                delete NowG;

            }

        } 

    }

    void ShadowFinder_DFS_Multithreading(tools::CGraph *  GlobalG,tools::Shadow * StartS,std::list<tools::Shadow *> & ShadowSave , int ThreadNum){

        std::vector<std::list<tools::Shadow *>> ThreadShadowSave(ThreadNum);

        std::list<std::thread> Threads;

        for(int _ =0 ; _ < ThreadNum ; _ ++){
            Threads.push_back(std::thread(ShadowFinder_DFS_Multithreading_Core,GlobalG,StartS,std::ref(ThreadShadowSave[_]),_,ThreadNum));
        }
        for(auto & _ : Threads){
            _.join();
        }
        for(int _ = 0 ; _ < ThreadNum ; _ ++){
            for(auto __ : ThreadShadowSave[_]){
                ShadowSave.push_back(__);
            }
        }

    }

}



namespace OnlineMT_Sample{

    bool OnceSampleMTSafe(tools::CGraph * OriginG , tools::Shadow * S){

        std::random_device rnd;
        std::mt19937 rd(rnd());
        auto a = new VertexIdx[S->nVerticesL+S->nVerticesR];
        // std::copy(S->VerticesL , S->VerticesL+S->nVerticesL , a); 
        std::memcpy(a,S->VerticesL,S->nVerticesL*sizeof(VertexIdx));
        // std::copy(S->VerticesR , S->VerticesR+S->nVerticesR , a+S->nVerticesL); 
        std::memcpy(a+S->nVerticesL, S->VerticesR, S->nVerticesR*sizeof(VertexIdx));
        std::shuffle(a, a+S->nVerticesL,rd); 
        std::shuffle(a+S->nVerticesL, a+S->nVerticesL+S->nVerticesR,rd); 

        auto re = IfConstructedClique(OriginG, a , a+S->nVerticesL, S->p, S->q);
        delete [] a ;
        return re;

    };

    double SampleShadow_Multithreading(tools::CGraph *  OriginG ,std::list<tools::Shadow *> & ShadowSave , double SWT , size_t  SampleTimes , size_t & NumThread,std::vector<double> & shadowWeight){
        
        std::vector<tools::Shadow *> ShadowVector;
        ShadowVector.reserve(ShadowSave.size());
        for(auto a : ShadowSave){
            ShadowVector.push_back(a);
        }

        auto __OnceSample_v2 = [&OriginG, &ShadowVector, &SampleTimes, &NumThread , &SWT, &shadowWeight ](double & totalWV , size_t ThreadId){
            sampleEngine onceIndex(shadowWeight);
            for(int _ = ThreadId ; _ < SampleTimes ; _+=NumThread){
                int index = onceIndex();
                if(OnceSampleMTSafe(OriginG, ShadowVector[index])){
                    totalWV+=SWT;
                }
            }
        }; 


        std::vector<double> ThreadtotalWV(NumThread,0);

        std::vector<std::thread> ThreadManage;

        for(int _ = 0 ; _ < NumThread ; _ ++) 
            ThreadManage.push_back(std::thread(__OnceSample_v2, std::ref(ThreadtotalWV[_]),_)); 
        for(auto & a : ThreadManage) 
            a.join(); 

        return std::accumulate(ThreadtotalWV.begin(),ThreadtotalWV.end(),(double)0.0);
        
    };

}



namespace local {

    tools::Shadow * SubGInduceByNodeR(tools::CGraph *OriginG,tools::Shadow * OriginS, VertexIdx Nodex){//

        if(Nodex < OriginG->nVerticesL) throw("check input");

        VertexIdx begin = 0 ;
        VertexIdx end = 0; 
        VertexIdx NumL =0 ;
        
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
        VertexIdx * & NewShadowVerticesL = TempNewShadowVerticesL;

        VertexIdx TempNumR = 0 ; 

        for(VertexIdx _  = 0 ; _ < NumL ; _ ++){
            OriginG->GetNebOfNode(NewShadowVerticesL[_], begin, end);
            TempNumR += end - begin;
        }


        VertexIdx * TempVerticesR = new VertexIdx[TempNumR];
        VertexIdx * it1 = TempVerticesR;

        for(VertexIdx _ = 0 ; _ < NumL ; _ ++){
            OriginG->GetNebOfNode(NewShadowVerticesL[_], begin, end);
            auto new_start  = tools::FindFirstGreaterThen(OriginG->col+begin, end - begin, Nodex);
            if(new_start == -1) continue;
            begin += new_start;
            for(VertexIdx __ = begin ; __ < end ; __ ++){
                *it1 = OriginG->col[__] ;
                it1++;
            }
        }

        if (TempVerticesR == it1) {
            delete[] TempNewShadowVerticesL;
            delete [] TempVerticesR;
            return NULL;

        }

        std::sort(TempVerticesR, it1);
        VertexIdx * it2 = TempVerticesR + 1 ; 
        for(VertexIdx * f = TempVerticesR + 1 ; f < it1 ; f++){
            if (*f == *(f-1)) {
                continue;
            } 
            else {
                *it2 = *f;
                it2++; 
            }
        }
        VertexIdx NumR = it2 - TempVerticesR;

        delete[] TempNewShadowVerticesL;
        delete [] TempVerticesR;

        auto aaaa =  new tools::Shadow({NumL,NumR,0,OriginS->p,OriginS->q-1,NULL,NULL});

        return aaaa;

    }

    void countP1InShadowSave(tools::CGraph * OriginG , std::list<tools::Shadow *> shadowSave,double & Z){


        auto iter = shadowSave.begin();

        while (iter!=shadowSave.end()) {
            if((*iter)->q == 1){
                for(int _  = 0 ; _ < (*iter)->nVerticesR ; _++){
                    if(OriginG->GetDegree((*iter)->VerticesR[_]) >= (*iter)->p){
                        Z+=boost::math::binomial_coefficient<double>(OriginG->GetDegree((*iter)->VerticesR[_]),(*iter)->p); 
                    }
                }
                iter = shadowSave.erase(iter);
            }
            else {
                iter++;
            }
        }

    }

    double doMission(tools::CGraph * OriginG, tools::Shadow * StartS,std::pair<int,int> Mission,size_t NumThread ){


        auto InduceByNodeR =  ShadowFinder(OriginG, StartS ,Mission.first);
        std::list<tools::Shadow *> ShadowSave;

        if (InduceByNodeR.tag == 0) {
            OnlineMT_Induce::ShadowFinder_DFS_Multithreading(OriginG, InduceByNodeR.G, ShadowSave, NumThread);
            delete InduceByNodeR.G;
        }
        else if(InduceByNodeR.tag == 1){
            ShadowSave.push_back(InduceByNodeR.G);
        }
        else {
            delete InduceByNodeR.G;
            return 0;
        }

        if (ShadowSave.size()==0) {
            return 0;
        }

        #ifdef MEMORY_USE

            size_t Temp ;
            auto b = [&ShadowSave,&Temp](){
                Temp = 0 ;
                for (auto & a : ShadowSave) {
                    Temp += (a->nVerticesL+a->nVerticesR)*4;
                }
                return Temp;
            };
            b();
            memCount::memoryTotal+=Temp;
            memCount::memoryMax = memCount::memoryMax > Temp ? memCount::memoryMax  : Temp;
            
        #endif
      


        double SWT = 0;

        std::vector<double> shadowWeights;
        shadowWeights.reserve(ShadowSave.size());
        for(auto eachShadow : ShadowSave){
            double SW = boost::math::binomial_coefficient<double>(eachShadow->nVerticesL, eachShadow->p) * boost::math::binomial_coefficient<double>(eachShadow->nVerticesR,eachShadow->q);
            SWT+=SW;
            shadowWeights.push_back(SW);
        } 

        auto & cumulative_distribution = shadowWeights;

        auto re = OnlineMT_Sample::SampleShadow_Multithreading(OriginG, ShadowSave, SWT,Mission.second, NumThread,cumulative_distribution);

        for(auto _  : ShadowSave){
            delete _;
        }

        return re;
    }

}



void OnlineSample_MT(tools::CGraph * OriginG,const int P , const int Q ,unsigned long long S , size_t NumThread){

    auto originGShadow = tools::csrToShadow(OriginG, P, Q);

    std::vector<double> phiV(originGShadow->nVerticesR,0);

    double phiT = 0;

    for(VertexIdx index = 0 ; index < originGShadow->nVerticesR ; index++){
        VertexIdx nodeRId = OriginG->nVerticesL + index ;
        phiV[index] = OriginG->GetDegree(nodeRId);
        phiT += phiV[index];
    }

    sampleEngine sampleBy_phiV(phiV);
    

    std::vector<unsigned int > nodeRSampleTime(OriginG->nVerticesR,0);
    for(unsigned long long _ = 0 ; _ < S ; _ ++){
        ++nodeRSampleTime[sampleBy_phiV()];
    }

    double W = 0;

    auto timeStart = std::chrono::system_clock::now();
    for(VertexIdx _ = 0; _ < OriginG->nVerticesR ;_++){

        if(nodeRSampleTime[_] == 0) continue;

        {
        #ifdef __DEBUG
            printf("all mission ：%d now mission ：%d  sample time:%u \n " , OriginG->nVerticesR , _  , nodeRSampleTime[_]);
            int ___ = _ + OriginG->nVerticesL;
        #endif
        }


    


        W+= (local::doMission(OriginG, originGShadow, {_+OriginG->nVerticesL,nodeRSampleTime[_]}, NumThread)/phiV[_]);
        
    }    

    {
        #ifdef  MEMORY_USE
            printf("all memory use  %lf max use %lf\n",(double)memCount::memoryTotal/(1024*1024) , (double)memCount::memoryMax /(1024*1024) );
            printf("all memory use %lu  max use  %lu\n",memCount::memoryTotal , memCount::memoryMax );
            printf("speed :%lf\n",((double)memCount::memoryTotal/(1024*1024)) / ((double)memCount::memoryMax /(1024*1024)));
        #endif
    }

    {
        #ifdef BOUND_COUNT
            size_t allShadowCount = boundCount::QEqual1 + boundCount::QEqual2 + boundCount::HitBound;
            printf("all shadow : %lu ,hit：%lu ,q=1 : %lu , q=2: %lu \n  hit rate ：%f ,q=1 rate: %f,q=2 rate: %f",
                    allShadowCount , boundCount::HitBound , boundCount::QEqual1 ,boundCount::QEqual2,
                    (double)boundCount::HitBound/allShadowCount,(double)boundCount::QEqual1/allShadowCount,(double)boundCount::QEqual2/allShadowCount);
            printf("\n");
            printf("\n");
        #endif
    }


    {   

        #ifdef SAMPLE_HIT
            extern   size_t hitTimes;
            extern   size_t allTime;
        #endif

        #ifdef SAMPLE_HIT
            printf("sample rate:%.5lf\n",(double)hitTimes/(double)(allTime));
            printf("sample times :%lu  , hit times:%lu\n",allTime,hitTimes);
            // printf("抽样命中率:%.5lf\n",(double)hitTimes/(double)(S));
        #endif

    }


    auto timeEnd = std::chrono::system_clock::now();

    printf("total: %.0lf\n" ,((W/(double)S)*phiT) );
    // std::cout<<"total: "<<(unsigned long long )((W/(double)S)*phiT)<<std::endl;
    auto allTime = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart);
    std::cout<<"time: " <<  (double)allTime.count() /1000 << std::endl;


}



