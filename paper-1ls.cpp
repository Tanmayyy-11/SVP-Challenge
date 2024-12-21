#include "paper-1ls.h"

template<class ZT, class FT>
Paper1_ls<ZT, FT>::Paper1_ls(const char *input_filename,int flags_bkz, int flags_gso, int prec, FloatType float_type)
: Paper1<ZT, FT>(input_filename,flags_bkz, flags_gso, prec, float_type) {
}


template<class ZT, class FT>
unique_ptr<Individual<ZT, FT>[]> Paper1_ls<ZT, FT>::runCrossMut(int dim, int k) {//local search
    unique_ptr<Individual<ZT,FT>[]> newPop = make_unique<Individual<ZT,FT>[]>(popObj->pop_size);
    newPop[0] = popObj->population[0];
    for(int i = 1;i<popObj->pop_size;i++)
    {
        unique_ptr<bool[]> a =  (popObj->encode(popObj->population[selection()].y.get(),popObj->totLength));
        unique_ptr<bool[]> b = (popObj->encode(popObj->population[selection()].y.get(),popObj->totLength));
        unique_ptr<bool[]> crossed = (cross(a.get(),b.get(),popObj->totLength));
        unique_ptr<bool[]> mutated  = (mutation(crossed.get(),popObj->totLength));
        newPop[i].y =  popObj->decode(mutated.get());
        newPop[i].x = newPop[i].YtoX(newPop[i].y.get(),popObj->get_mu(),popObj->dim);
        newPop[i].norm = newPop[i].get_norm(newPop[i].matrix_multiply(newPop[i].x.get(),popObj->get_B(),popObj->dim).get(),popObj->dim);
        newPop[i].dim = dim;
        if(newPop[i].norm==0)
        {
            i--;
            continue;
        }
        int p=1;//p=1 => last iteration of local search had an improvement
        
        while(p==1){
            int pos=-1;
            FT norm =newPop[i].norm;
            p=0;
            int sgn;
            for(int j=0;j<dim;j++)
            {
                ZT diff ;
                diff.sub_ui(newPop[i].y[j],1);
                diff.abs(diff);
                FT absdiff;
                absdiff.set_z(diff);
                if(absdiff <= popObj->alpha[j]){
                    newPop[i].y[j].sub_ui(newPop[i].y[j],1);
                    unique_ptr<ZT[]> use=newPop[i].YtoX(newPop[i].y.get(),popObj->get_mu(),popObj->dim);
                    unique_ptr<ZT[]> result = newPop[i].matrix_multiply(use.get(),popObj->get_B(),popObj->dim);
                    FT temp=newPop[i].get_norm(result.get(),popObj->dim);
                    if((temp!=0)&&(temp<norm))
                    {
                        pos=j;
                        sgn=-1;
                        norm=temp;
                    }
                    newPop[i].y[j].add_ui(newPop[i].y[j],1);
                }
                ZT sum ;
                sum.add_ui(newPop[i].y[j],1);
                sum.abs(sum);
                FT abssum;
                abssum.set_z(sum);
                if(abssum <= popObj->alpha[j]){
                    newPop[i].y[j].add_ui(newPop[i].y[j],1);
                    unique_ptr<ZT[]> use2=newPop[i].YtoX(newPop[i].y.get(),popObj->get_mu(),popObj->dim);
                    unique_ptr<ZT[]> result2 = newPop[i].matrix_multiply(use2.get(),popObj->get_B(),popObj->dim);
                    FT temp =newPop[i].get_norm(result2.get(),popObj->dim);
                    if((temp!=0)&&(temp<norm))
                    {
                        pos=j;
                        sgn=+1;
                        norm=temp;
                    }
                    newPop[i].y[j].sub_ui(newPop[i].y[j],1);
                }
            }
            if(pos!=(-1))
            {
                // newPop[i].y[pos]+=sgn;
                p=1;//improvement occured
                // cout<<"local search worked"<<endl;
                if(sgn==-1)
                {
                    newPop[i].y[pos].sub_ui(newPop[i].y[pos],1);
                }
                else if(sgn==1)
                {
                    newPop[i].y[pos].add_ui(newPop[i].y[pos],1);

                }
                newPop[i].x=newPop[i].YtoX(newPop[i].y.get(),popObj->get_mu(),popObj->dim);
                newPop[i].norm=norm;
            }
        }
        
    }
    return newPop;
}

template class Paper1_ls<Z_NR<mpz_t>,FP_NR<mpfr_t>>;