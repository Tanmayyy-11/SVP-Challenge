#include "paper-1.h"
using namespace std;
using namespace fplll;
template<class ZT, class FT>
class Paper1_ls:public Paper1<ZT,FT>
{
    public:
    using GA<ZT, FT>::popObj;
    using Paper1<ZT, FT>::selection;
    using Paper1<ZT, FT>::cross;
    using Paper1<ZT, FT>::mutation;
    using Paper1<ZT, FT>::runCrossMut;
    using Paper1<ZT, FT>::compare;
    using Paper1<ZT, FT>::logical_xor;
    using Paper1<ZT, FT>::logical_and;
    using Paper1<ZT, FT>::runGA;
    Paper1_ls(const char *input_filename,int flags_bkz, int flags_gso, int prec, FloatType float_type);
    // Paper1_ls<ZT, FT>::Paper1_ls(const char *input_filename, int flags_bkz, int flags_gso, int prec, FloatType float_type): Paper1<ZT, FT>(input_filename, flags_bkz, flags_gso, prec, float_type)
    unique_ptr<Individual<ZT, FT>[]> runCrossMut(int dim, int k) override ;
};