#include <algorithm>
#include <ios>
#include <iostream>
#include "graph.h"
#include "sample.h"
#include "onlinesample.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <istream>
#include <iterator>
#include <list>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>
#include <random>
#include "ultis.hpp"
#include <cstdio>
#define DEBUG

using namespace std;

        


int main(int argc , char * argv[]){

    if (argc < 2) {
        cout<<"check input\n";
        exit(0);
    }

    {//online sample
        SpecialBigraph a(argv[1],stoi(argv[2]),stoi(argv[3]),0) ;
        printf(" ------  P = %d  Q = %d \n",stoi(argv[2]),stoi(argv[3]));
        a.read_graph();
        a.prepare_graph();  
        auto b = a.sbgraphToShadow();
        OnlineSample_MT(b, b->p, b->q, stoi(argv[4]),stoi(argv[5]));
        delete b; 
    }




}
