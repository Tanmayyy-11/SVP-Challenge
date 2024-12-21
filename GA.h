#pragma once
#include "repres.h"
using namespace std;
using namespace fplll;
template<class ZT,class FT> class GA
{
    public: 
    unique_ptr<repres<ZT,FT>> popObj;
    // virtual int selection(int k)=0;
    virtual int selection()=0;
    virtual unique_ptr<bool[]>  cross(bool* a,bool* b,ZT tot_length)=0;
    virtual unique_ptr<bool[]> mutation(bool* a,ZT tot_length)=0;
    GA(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type);
    virtual unique_ptr<ZT[]>  runGA(FT targetNorm,int k,Individual<ZT,FT>vb)=0;
    virtual unique_ptr<Individual<ZT,FT>[]>runCrossMut(int dim,int k)=0;
};
