#pragma once
#include <iostream>
#include <fplll.h>
#include <cmath>
#include <fplll.h>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include<random>
#include "test_utils.h"
using namespace std;
using namespace fplll;
#include "GA.h"

template<class ZT, class FT>
class Paper1: public GA<ZT, FT> {
public:
    using GA<ZT, FT>::popObj;
    // int selection(int k) override;
    int selection()override;
    unique_ptr<bool[]> cross(bool* a, bool* b, ZT tot_length) override;
    unique_ptr<bool[]> mutation(bool* a, ZT tot_length) override;
    unique_ptr<Individual<ZT,FT>[]>  runCrossMut(int dim, int k) override;
    static bool compare(Individual<ZT, FT> &i1, Individual<ZT, FT> &i2);
    Paper1(const char *input_filename,int flags_bkz, int flags_gso, int prec, FloatType float_type);
    unique_ptr<bool[]> logical_xor(bool* a, bool* b, ZT tot_length);
    unique_ptr<bool[]> logical_and(bool* a, bool* b, ZT tot_length);
    unique_ptr<ZT[]>  runGA(FT targetNorm, int k,Individual<ZT,FT>vb) override;
};
