#pragma once
#include <cmath>
#include <fplll.h>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include<random>
#include "test_utils.h"
#include "preprocessing.h"
#include "Individual.h"
using namespace std;
using namespace fplll;

template<class ZT,class FT>
class repres {
public:
    int pop_size;
    int dim;
    unique_ptr<FT[]> alpha;
    unique_ptr<ZT[]> length;
    ZT totLength;
    unique_ptr<Individual<ZT,FT>[]> population;
    unique_ptr<preProcessing<ZT,FT>> preprocess;
    unique_ptr<ZT[]>  decode(bool *chromosome);
    FT get_norm(ZT* vect, int dim);
    FT get_norm(FT* vect, int dim);
    unique_ptr<bool[]>  encode(ZT* y, ZT totalLength);
    void initialise(Individual<ZT,FT>v0);
    repres(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type);
    unique_ptr<ZT[]>* get_B();
    unique_ptr<FT[]>* get_mu();
};
