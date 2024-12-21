#include"paper-1.h"
using namespace std;
using namespace fplll;
// template<class ZT, class FT>
// int Paper1<ZT, FT>::selection(int k) {
//     int pop_size=popObj->pop_size;
//     unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//     std::mt19937 generator(seed);
//     std::uniform_int_distribution<int> distribution(0, pop_size - 1);
//     FT comp;
//     int val=0;
//     for(int i=0;i<k;i++)
//     {
//         int ind=distribution(generator);
//         FT use=(popObj->population[ind]).norm;
//         if((comp>use)||(i==0))
//         {
//             comp=use;
//             val=ind;
//         }
//     }
//     return val;
// }
template<class ZT, class FT>
int Paper1<ZT, FT>::selection() {
    int pop_size=popObj->pop_size;
    FT max_norm;
    for(int i=0;i<pop_size;i++)
    {
        FT curr_norm=(popObj->population[i]).norm;
        if((i==0)||(max_norm>curr_norm))max_norm=curr_norm;
    }
    max_norm=(max_norm*max_norm)*max_norm;
    // auto fitness_values{make_unique<ZT[]>(pop_size)};
    ZT totalsum;
    totalsum=0;
    for(int i=0;i<pop_size;i++)
    {
        FT scaled_norm=max_norm/((popObj->population[i]).norm);
        scaled_norm.floor(scaled_norm);
        ZT value;
        value.set_f(scaled_norm);
        totalsum.add(totalsum,value);
    }
    ZT randomNumber;
    randomNumber.randm(totalsum);
    int index=pop_size-1;
    for(int i=0;i<pop_size;i++)
    {
        FT scaled_norm=max_norm/((popObj->population[i]).norm);
        scaled_norm.floor(scaled_norm);
        ZT value;
        value.set_f(scaled_norm);
        randomNumber.sub(randomNumber,value);
        ZT comp;
        comp=0;
        if(randomNumber<comp)return i;
    }
    return index;
}
template<class ZT, class FT>
unique_ptr<bool[]> Paper1<ZT, FT>::cross(bool* a, bool* b, ZT tot_length) {
     random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distribution(0, 1);
    unique_ptr<bool[]>randomu=make_unique<bool[]>(int(tot_length.get_d()));
    unique_ptr<bool[]>randomubar=make_unique<bool[]>(int(tot_length.get_d()));
    for(int i=0;i<int(tot_length.get_d());i++)
    {
        int random_val=distribution(gen);
        if(random_val)
        {
            randomu[i]=false;
            randomubar[i]=true;
        }
        else{
            randomu[i]=true;
            randomubar[i]=false;
        }
    }
    unique_ptr<bool[]> ans=  logical_xor(logical_and(a,randomu.get(),tot_length).get(),logical_and(b,randomubar.get(),tot_length).get(),tot_length);
    // delete[]randomu;
    // delete[] randomubar;
    return ans;
}

template<class ZT, class FT>
unique_ptr<bool[]> Paper1<ZT, FT>::mutation(bool* a, ZT tot_length) {
    unique_ptr<bool[]> randomM=make_unique<bool[]>(int(tot_length.get_d()));
   random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distribution(0, ((int)tot_length.get_d())-1);
   for(int i=0;i<int(tot_length.get_d());i++)
   {
    int random_val=distribution(gen);
    if(random_val<1)randomM[i]=true;
    else randomM[i]=false;
   }
   unique_ptr<bool[]>ans= logical_xor(a,randomM.get(),tot_length);
//    delete[] randomM;
   return ans;
}
template<class ZT, class FT>
unique_ptr<bool[]> Paper1<ZT, FT>::logical_xor(bool* a, bool* b, ZT tot_length) {
     unique_ptr<bool[]> result = make_unique<bool[]>(int(tot_length.get_d()));
    for (int i = 0; i < int(tot_length.get_d()); ++i) {
        result[i] = a[i] ^ b[i];
    }
    return result;
}

template<class ZT, class FT>
unique_ptr<bool[]> Paper1<ZT, FT>::logical_and(bool* a, bool* b, ZT tot_length) {
     unique_ptr<bool[]> result = make_unique<bool[]>(int(tot_length.get_d()));
    for (int i = 0; i < int(tot_length.get_d()); ++i) {
        result[i] = a[i] && b[i];
    }
    return result;
}
template<class ZT, class FT>
unique_ptr<Individual<ZT,FT>[]> Paper1<ZT, FT>::runCrossMut(int dim, int k) {
     unique_ptr<Individual<ZT,FT>[]> newPop = make_unique<Individual<ZT,FT>[]>(popObj->pop_size);
    // newPop[0] = new Individual<ZT,FT>(dim);
    // for(int i = 0;i<dim;i++)
    // {
    //     newPop[0].x[i] = popObj->population[0].x[i];
    //     newPop[0].y[i] = popObj->population[0].y[i];
    // }
    // newPop[0].norm = popObj->population[0].norm;
    // for(int i = 0;i<dim;i++)
    // {
    //     newPop[0].x[i] = popObj->population[0].x[i];
    //     newPop[0].y[i] = popObj->population[0].y[i];
    // }
    // newPop[0].norm = popObj->population[0].norm;
    newPop[0] = popObj->population[0];
    // cout<<"Herisfine\n";
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
            // delete[] newPop[i].y;
            // delete[] newPop[i].x;
            i--;
        }
       // cout<<i<<" "<<newPop[i].norm<<'\n';
        // delete[]a;
        // delete[]b;
        // delete[]crossed;
        // delete[]mutated;
    }
    return newPop;
}
template<class ZT, class FT>
bool Paper1<ZT, FT>::compare(Individual<ZT, FT>& i1, Individual<ZT, FT>& i2) {
    return (i1.norm<i2.norm);
}

template<class ZT, class FT>
Paper1<ZT, FT>::Paper1(const char *input_filename,int flags_bkz, int flags_gso, int prec, FloatType float_type)
    : GA<ZT, FT>(input_filename,flags_bkz, flags_gso, prec, float_type) {
}

template<class ZT, class FT>
unique_ptr<ZT[]>  Paper1<ZT, FT>::runGA(FT targetNorm, int k,Individual<ZT,FT> vb) {
    popObj->initialise(vb);
    // cout<<"Initialisation Complete\n";
    sort(popObj->population.get(),popObj->population.get()+popObj->pop_size,compare);
    // cout<<popObj->population[0].norm<<" "<<popObj->population[1].norm<<" "<<popObj->population[2].norm<<"\n";
    // Individual<ZT,FT>v0;
    // for(int i = 0;i<popObj->dim;i++)
    // {
    //     v0.x[i] = popObj->population[0].x[i];
    //     v0.y[i] = popObj->population[0].y[i];
    // }
    // v0.norm = popObj->population[0].norm;
    // cout<<"Sorting Complete\n";
    Individual<ZT,FT> v0 (popObj->population[0]);
    // cout<<"Copy Working\n";
    int iter=0;
    FT prevNorm;
    prevNorm = 0.0;
    cout<<"norm v0 "<<v0.norm<<endl;
    cout<<"target norm "<< targetNorm <<endl;
    ZT a,b,c,d;
    a=-1;b=0;c=1;d=0;
    ZT dummy_x[4]={a,b,c,d};
    unique_ptr<ZT[]>yy=v0.XtoY(dummy_x,popObj->preprocess->mu.get(),4);//debug
    cout<<"OUR Y"<<endl;
    for(int i=0;i<4;i++)cout<<yy[i]<<" ";
    cout<<endl;
    while(v0.norm>targetNorm)
    {
        // cout<<"Inside For loop\n";
        popObj->population = runCrossMut(popObj->dim,k);
        // cout<<"CrossMut done\n";
        sort(popObj->population.get(),popObj->population.get()+popObj->pop_size,compare);
        if((iter%1000)==0)//debug
        {
        cout<<iter<<"\n";
        cout<<v0.norm<<" "<<prevNorm<<"\n";
        }
    //    for(int i = 0;i<popObj->dim;i++)
    //     {
    //         v0.x[i] = popObj->population[0].x[i];
    //         v0.y[i] = popObj->population[0].y[i];
    //     }
    //     v0.norm = popObj->population[0].norm;
    //     v0.dim = popObj->po
        v0 = popObj->population[0];

        if(prevNorm!=v0.norm)
        {
             cout<<v0.norm<<" "<<'\n';
             prevNorm = v0.norm;
             iter = 0;
        }
        // iter.add_ui(iter,1);
        iter++;//debug
        // if(iter==100000)//debug
        // {
        //     cout<<"No. of Iteration exceeds 1000\n";
        //     cout<<"Restarting Algorithm\n";
        //     // delete []popObj->population;
        //     return runGA(targetNorm,k,v0);
        // }

    }
    // cout<<1<<endl;
    // for(e:v0.y)cout<<e<<" ";
    // cout<<endl;
    for(int i=0;i<4;i++)
    {
        cout<<v0.x.get()[i]<<" ";
    }
    cout<<endl;
    return v0.matrix_multiply(v0.x.get(),popObj->get_B(),popObj->dim);
}

template class Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>;
