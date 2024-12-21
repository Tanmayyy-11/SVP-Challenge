#pragma once
#include <cmath>
#include <fplll.h>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <memory>
#include<random>
#include "test_utils.h"
using namespace std;
using namespace fplll;
template<class ZT,class FT>class preProcessing
{
    public:
        int dim;
        Matrix<ZT>U;
        Matrix<ZT>UT;
        std::unique_ptr<std::unique_ptr<ZT[]>[]> B;
        std::unique_ptr<std::unique_ptr<FT[]>[]> Bstar;
        std::unique_ptr<std::unique_ptr<FT[]>[]> mu;
        preProcessing(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type);
        preProcessing(unique_ptr<std::unique_ptr<ZT[]>[]>& old_B,unique_ptr<ZT[]>& vect,int dimension,int flags_gso);
};