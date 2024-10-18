#ifndef ONLINESAMPLE
#define ONLINESAMPLE
#include "graph.h"

void OnlineSample(tools::CGraph * OriginG,const int P , const int Q ,unsigned long long S);

void OnlineSample_MT(tools::CGraph * OriginG,const int P , const int Q ,unsigned long long S , size_t NumThread);




#endif






