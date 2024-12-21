// #pragma once
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
using namespace std::chrono;

template<typename T>
void writeMatrixToFile(ZZ_mat<T>& matrix, const std::string& filename) {
    std::ofstream outFile(filename);  // Open file for writing

    if (!outFile) {  // Check if the file opened successfully
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    int dim = matrix.get_rows();  // Get the number of rows
    int cols = matrix.get_cols();
    outFile << dim << " " << cols << std::endl;  // Write dimensions

    for (int i=0;i<dim;i++) {  // Iterate through each row
        for (int j=0;j<cols;j++) {  // Iterate through each element in the row
            outFile << matrix[i][j] << " ";  // Write the element
        }
        outFile << std::endl;  // New line after each row
    }

    outFile.close();  // Close the file
    std::cout << "Matrix written to " << filename << std::endl;
}


template<class ZT,class FT>class preProcessing
{
    public:
        int dim;
        Matrix<ZT>U;
        Matrix<ZT>UT;
        ZT** B;
        FT** Bstar;
        FT** mu;
        preProcessing(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type);
};
template<class ZT,class FT>preProcessing<ZT,FT>::preProcessing(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type)
{
    ZZ_mat<mpz_t>A;
    int status = read_file(A, input_filename);
    dim = A.get_rows();
    //cout<<dim<<endl;
    if(status==0)
    {
        int status = bkz_reduction(A, 8, flags_bkz, float_type,prec);
        MatGSO<ZT, FT> M(A, U, UT, flags_gso);
        M.update_gso();
        writeMatrixToFile(A,"inputbases.txt");
        B = new ZT*[dim];
        Bstar = new FT*[dim];
        mu = new FT*[dim];
        // alpha = new FT[dim];
        // length = new ZT[dim];
        for(int i = 0;i<dim;i++)
        {
            Bstar[i] = new FT[dim];
            mu[i] = new FT[dim];
            B[i] = new ZT[dim];
        }
        for(int i = 0;i<dim;i++)
        {
            for(int j = 0;j<dim;j++)
            {
                M.get_gram(Bstar[i][j],i,j);
                M.get_mu(mu[i][j],i,j);
                // cout<<mu[i][j]<<'\n';
                B[i][j] = A[i][j];
            }
        }
        FT norm = 0.0;
        for (int i = 0; i < dim; ++i) {
            // FT use=B[0][i];
            FT use;
            use.set_z(B[0][i]);
            norm += use * use;
        }
        cout<<"first vec norm "<< sqrt(norm)<<endl;
        // cout<<B[i].get
        // cout<<1<<endl;
    }
    else
    {
        cout<<"Error in read file\n";
        
    }


}
template<class ZT,class FT>
class Individual {
public:
    FT norm;
    static  ZT* common_memory; // Initialization with nullptr
    ZT* x;
    ZT* y;
    FT get_norm(ZT* vect, int dim);
    ZT* matrix_multiply(ZT* vec, ZT** mat, int dim);
    ZT* YtoX(ZT* y,FT**mu,int dim);
    Individual(int dim, FT** mu, FT* alpha, ZT** B, FT** Bstar);
    Individual();
    ~Individual();
};

// // Definition of static member outside the class
template<class ZT, class FT>
ZT* Individual<ZT, FT>::common_memory = new ZT[1];
template<class ZT,class FT>
ZT* Individual<ZT,FT>::matrix_multiply(ZT* vec, ZT** mat, int dim) {
    ZT* result = new ZT[dim];
    for (int i = 0; i < dim; ++i) {
        result[i] = 0;
        for (int j = 0; j < dim; ++j) {
            ZT temp; 
            temp.mul(vec[j], mat[i][j]);
            result[i].add(result[i], temp);
        }
    }
    return result;
}
template<class ZT,class FT>
ZT* Individual<ZT,FT>::YtoX(ZT* y,FT** mu, int dim) {
    ZT* x = new ZT[dim];
    for (int i = dim - 1; i >= 0; i--) {
        FT t=0.0;
        for(int j=i+1;j<dim;j++)
        {
            FT temp;
            temp.set_z(x[j]);
            t+=mu[j][i]*temp;
        }
        t.rnd(t);
        ZT t_z;
        t_z.set_f(t);
        x[i].sub(y[i],t_z);
    }
    return x;
}
template<class ZT,class FT>
Individual<ZT,FT>::Individual(int dim, FT** mu, FT* alpha, ZT** B, FT** Bstar) {

    this->x = new ZT[dim];
    this->y = new ZT[dim];
    ZT* alpha_r=new ZT[dim];
    for(int i=0;i<dim;i++)
    {
        FT temp=alpha[i];
        alpha[i].floor(alpha[i]);
        ZT use;
        use.set_f(alpha[i]);
        alpha_r[i].add_ui(use,1);
        alpha[i]=temp;
    }
    for (int i = dim - 1; i >= 0; i--) {
        y[i].randm(alpha_r[i]);
        cout<<y[i] <<" ";
        FT t=0.0;
        for(int j=i+1;j<dim;j++)
        {
            FT temp;
            temp.set_z(x[j]);
            t+=mu[j][i]*temp;
        }
        t.rnd(t);
        ZT t_z;
        t_z.set_f(t);
        x[i].sub(y[i],t_z);    }
        cout<<endl;
    ZT* vect = matrix_multiply(x, B, dim);
    this->norm = get_norm(vect, dim);
    delete[] vect;
}
template<class ZT,class FT>
Individual<ZT,FT>::~Individual() 
{
}
template<class ZT,class FT>
Individual<ZT,FT>::Individual()
{
    this->x = NULL;
    this->y = NULL;
    this->norm = 0.0;
    // this->x=common_memory;
    // this->y=common_memory;
    // x=new ZT[]
}
template<class ZT,class FT>
class repres {
public:
    int pop_size;
    int dim;
    FT* alpha;
    ZT* length;
    ZT totLength;
    Individual<ZT,FT>* population;
    preProcessing<ZT,FT>* preprocess;
    ZT* decode(bool *chromosome);
    FT get_norm(ZT* vect, int dim);
    FT get_norm(FT* vect, int dim);
    bool* encode(ZT* y, ZT totalLength);
    // void initialise();
    void initialise(Individual<ZT,FT>v0);
    repres(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type);
    ZT **get_B();
    FT** get_mu();
};
template<class ZT,class FT>repres<ZT,FT>::repres(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type)
{

    preprocess = new preProcessing<ZT,FT>(input_filename, flags_bkz, flags_gso, prec, float_type);
    dim = preprocess->dim;
    pop_size = 2*dim;
    FT normB0 = get_norm(preprocess->Bstar[0],dim);
    FT log2 = log(2);
    alpha=new FT[dim];
    length=new ZT[dim];
    totLength = 0;
    for(int i = 0;i<dim;i++)
    {
        alpha[i] = normB0/get_norm(preprocess->Bstar[i],dim);
        length[i] .set_f(floor(log(alpha[i])/log2)+FT(2));
        totLength.add(totLength,length[i]);
    }


}
template <class ZT,class FT>
FT Individual<ZT,FT>::get_norm(ZT* vect, int dim) {
    FT norm = 0.0;
    for (int i = 0; i < dim; ++i) {
        FT use;
        use.set_z(vect[i]);
        norm += use * use;
    }
    return sqrt(norm);
}
template <class ZT,class FT>
FT repres<ZT,FT>::get_norm(FT* vect, int dim) {
    FT norm = 0.0;
    for (int i = 0; i < dim; ++i) {
        FT use=vect[i];
        norm += use * use;
    }
    return sqrt(norm);
}

template <class ZT,class FT>
FT repres<ZT,FT>::get_norm(ZT* vect, int dim) {
    FT norm = 0.0;
    for (int i = 0; i < dim; ++i) {
        FT use = vect[i].get_d();
        norm += use * use;
    }
    return sqrt(norm);
}
template<class ZT,class FT>
ZT* repres<ZT,FT>::decode(bool *chromosome) {
    ZT *y = new ZT[dim];
    int start = 0;
    for (int i = 0; i < dim; i++) {
        ZT mult;
        mult=(long)1; // Changed pointer to value
        for (int j = start + (int)length[i].get_d() - 1; j > start; j--) {
            ZT use = mult;
            long temp=chromosome[j];
            use.mul_si(use,temp);
            y[i].add(y[i], use);
            mult.mul_si(mult, (long)2);
        }
        if (chromosome[start]) y[i].mul_si(y[i], (long)-1);
        start = start + (int)length[i].get_d();
    }
    return y;
}

template<class ZT,class FT>
bool* repres<ZT,FT>::encode(ZT* y, ZT totalLength) {
    bool* chromosome = new bool[(int)totalLength.get_d()];
    ZT tot;
    tot = (long)0;
    for (int i = 0; i < dim; i++) {
        bool sign = 0;
        ZT element = y[i];
        if (y[i] < (long)0) {
            sign = 1;
            element.neg(element);
        }
        ZT curEl = y[i];
        for (int j = tot.get_d() + length[i].get_d() - 1; j > tot.get_d(); j--) {
            ZT curBit, mod_;
            mod_=(long)2;
            curBit.mod(curEl, mod_);
            chromosome[j] = (int)curBit.get_d();
            curEl.div_2si(curEl, 1);//divide by 2^1
        }
        chromosome[(int)tot.get_d()] = sign;
        tot.add(tot, length[i]);
    }
    return chromosome;
}

template<class ZT,class FT>
void repres<ZT,FT>::initialise(Individual<ZT,FT>v0) {
    //cout<<pop_size<<endl;

    population = new Individual<ZT,FT>[pop_size];
    if(v0.x==NULL)
    {
        population[0].y = new ZT[dim] ;
        population[0].y[0] = 1;
        for(int i = 1;i<dim;i++)population[0].y[i] = 0;

        population[0].x = population[0].YtoX(population[0].y,preprocess->mu,dim);
        ZT* vect = population[0].matrix_multiply(population[0].x, preprocess->B,dim);
        population[0].norm = population[0].get_norm(vect,dim);
    }
    else 
        population[0] = v0;
    for (int i = 1; i < pop_size; i++) {
       // cout<<1<<endl;
        population[i] = Individual<ZT,FT>(dim, preprocess->mu, alpha, preprocess->B, preprocess->Bstar);
    }
}
template<class ZT,class FT>ZT** repres<ZT,FT>::get_B()
{
    return preprocess->B;
}
template<class ZT,class FT>FT** repres<ZT,FT>::get_mu()
{
    return preprocess->mu;
}
template<class ZT,class FT> class GA
{
    public: 
    repres<ZT,FT> *popObj;
    virtual int selection(int k)=0;
    virtual bool* cross(bool* a,bool* b,ZT tot_length)=0;
    virtual bool* mutation(bool* a,ZT tot_length)=0;
    GA(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type);
    virtual ZT* runGA(FT targetNorm,int k,Individual<ZT,FT>vb)=0;
    virtual Individual<ZT,FT>*runCrossMut(int dim,int k)=0;
};
template<class ZT, class FT>
GA<ZT, FT>::GA(const char *input_filename,int flags_bkz, int flags_gso, int prec, FloatType float_type) {
    popObj = new repres<ZT, FT>(input_filename,flags_bkz, flags_gso, prec, float_type);
}
template<class ZT, class FT>
class Paper1: public GA<ZT, FT> {
public:
    using GA<ZT, FT>::popObj;
    int selection(int k) override;
    bool* cross(bool* a, bool* b, ZT tot_length) override;
    bool* mutation(bool* a, ZT tot_length) override;
    Individual<ZT, FT>* runCrossMut(int dim, int k) override;
    static bool compare(Individual<ZT, FT> &i1, Individual<ZT, FT> &i2);
    Paper1(const char *input_filename,int flags_bkz, int flags_gso, int prec, FloatType float_type);
    bool* logical_xor(bool* a, bool* b, ZT tot_length);
    bool* logical_and(bool* a, bool* b, ZT tot_length);
    ZT* runGA(FT targetNorm, int k,Individual<ZT,FT>vb) override;
};
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
    Individual<ZT, FT>* runCrossMut(int dim, int k) override;
};
template<class ZT, class FT>
int Paper1<ZT, FT>::selection(int k) {
    int pop_size=popObj->pop_size;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_int_distribution<int> distribution(0, pop_size - 1);
    FT comp;
    int val=0;
    for(int i=0;i<k;i++)
    {
        int ind=distribution(generator);
        FT use=(popObj->population[ind]).norm;
        if((comp>use)||(i==0))
        {
            comp=use;
            val=ind;
        }
    }
    return val;
}
template<class ZT, class FT>
bool* Paper1<ZT, FT>::cross(bool* a, bool* b, ZT tot_length) {
     random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distribution(0, 1);
    bool *randomu=new bool[int(tot_length.get_d())];
    bool *randomubar=new bool[int(tot_length.get_d())];
    for(int i=0;i<int(tot_length.get_d());i++)
    {
        int random_val=distribution(gen);
        if(random_val)
        {
            randomu[i]=false;
            randomubar[i]=true;
        }
        else{
            randomu[i]=false;
            randomubar[i]=true;
        }
    }
    bool* ans=  logical_xor(logical_and(a,randomu,tot_length),logical_and(b,randomubar,tot_length),tot_length);
    delete[]randomu;
    delete[] randomubar;
    return ans;
}

template<class ZT, class FT>
bool* Paper1<ZT, FT>::mutation(bool* a, ZT tot_length) {
    bool* randomM=new bool[int(tot_length.get_d())];
   random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> distribution(0, ((int)tot_length.get_d())-1);
   for(int i=0;i<int(tot_length.get_d());i++)
   {
    int random_val=distribution(gen);
    if(random_val<1)randomM[i]=true;
    else randomM[i]=false;
   }
   bool*ans= logical_xor(a,randomM,tot_length);
   delete[] randomM;
   return ans;
}
template<class ZT, class FT>
bool* Paper1<ZT, FT>::logical_xor(bool* a, bool* b, ZT tot_length) {
     bool* result = new bool[int(tot_length.get_d())];
    for (int i = 0; i < int(tot_length.get_d()); ++i) {
        result[i] = a[i] ^ b[i];
    }
    return result;
}

template<class ZT, class FT>
bool* Paper1<ZT, FT>::logical_and(bool* a, bool* b, ZT tot_length) {
    bool* result = new bool[int(tot_length.get_d())];
    for (int i = 0; i < int(tot_length.get_d()); ++i) {
        result[i] = a[i] && b[i];
    }
    return result;
}
template<class ZT, class FT>
Individual<ZT, FT>* Paper1<ZT, FT>::runCrossMut(int dim, int k) {
     Individual<ZT,FT>* newPop = new Individual<ZT,FT>[popObj->pop_size];
    newPop[0] = popObj->population[0];
    //cout<<popObj->pop_size<<'\n';
    for(int i = 1;i<popObj->pop_size;i++)
    { 
        bool *a = popObj->encode(popObj->population[selection()].y,popObj->totLength);
        bool *b = popObj->encode(popObj->population[selection()].y,popObj->totLength);
        bool *crossed = cross(a,b,popObj->totLength);
        bool *mutated = mutation(crossed,popObj->totLength);
        newPop[i].y =  popObj->decode(mutated);
        newPop[i].x = newPop[i].YtoX(newPop[i].y,popObj->get_mu(),popObj->dim);
        newPop[i].norm = newPop[i].get_norm(newPop[i].matrix_multiply(newPop[i].x,popObj->get_B(),popObj->dim),popObj->dim);
        if(newPop[i].norm==0)
        {
            delete[] newPop[i].y;
            delete[] newPop[i].x;
            i--;
        }
       // cout<<i<<" "<<newPop[i].norm<<'\n';
        delete[]a;
        delete[]b;
        delete[]crossed;
        delete[]mutated;
    }
    return newPop;
}
template<class ZT, class FT>
Individual<ZT, FT>* Paper1_ls<ZT, FT>::runCrossMut(int dim, int k) {
     Individual<ZT,FT>* newPop = new Individual<ZT,FT>[popObj->pop_size];
    newPop[0] = popObj->population[0];
    for(int i = 1;i<popObj->pop_size;i++)
    {
        bool *a = popObj->encode(popObj->population[selection()].y,popObj->totLength);
        bool *b = popObj->encode(popObj->population[selection()].y,popObj->totLength);
        bool *crossed = cross(a,b,popObj->totLength);
        bool *mutated = mutation(crossed,popObj->totLength);
        newPop[i].y =  popObj->decode(mutated);
        newPop[i].x = newPop[i].YtoX(newPop[i].y,popObj->get_mu(),popObj->dim);
        newPop[i].norm = newPop[i].get_norm(newPop[i].matrix_multiply(newPop[i].x,popObj->get_B(),popObj->dim),popObj->dim);
        if(newPop[i].norm==0)
        {
        delete[] newPop[i].y;
        delete[] newPop[i].x;
        i--;
        delete[]a;
        delete[]b;
        delete[]crossed;
        delete[]mutated;
        continue;
        }
        FT norm =newPop[i].norm;
        int pos=-1;
        int sgn;
        for(int j=0;j<dim;j++)
        {
            newPop[i].y[j].sub_ui(newPop[i].y[j],1);
            ZT* use=newPop[i].YtoX(newPop[i].y,popObj->get_mu(),popObj->dim);
            FT temp=newPop[i].get_norm(newPop[i].matrix_multiply(use,popObj->get_B(),popObj->dim),popObj->dim);
            if((temp!=0)&&(temp<norm))
            {
                pos=j;
                sgn=-1;
                norm=temp;
            }
            newPop[i].y[j].add_ui(newPop[i].y[j],2);
            use=newPop[i].YtoX(newPop[i].y,popObj->get_mu(),popObj->dim);
            temp=newPop[i].get_norm(newPop[i].matrix_multiply(use,popObj->get_B(),popObj->dim),popObj->dim);
            if((temp!=0)&&(temp<norm))
            {
                pos=j;
                sgn=+1;
                norm=temp;
            }
            newPop[i].y[j].sub_ui(newPop[i].y[j],1);
        }
        if(pos!=(-1))
        {
            // newPop[i].y[pos]+=sgn;
            if(sgn==-1)
            {
                newPop[i].y[pos].sub_ui(newPop[i].y[pos],1);
            }
            else if(sgn==1)
            {
                newPop[i].y[pos].add_ui(newPop[i].y[pos],1);

            }
            newPop[i].x=newPop[i].YtoX(newPop[i].y,popObj->get_mu(),popObj->dim);
            newPop[i].norm=norm;
        }
        delete[]a;
        delete[]b;
        delete[]crossed;
        delete[]mutated;
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
Paper1_ls<ZT, FT>::Paper1_ls(const char *input_filename,int flags_bkz, int flags_gso, int prec, FloatType float_type)
    : Paper1<ZT, FT>(input_filename,flags_bkz, flags_gso, prec, float_type) {
}

template<class ZT, class FT>
ZT* Paper1<ZT, FT>::runGA(FT targetNorm, int k,Individual<ZT,FT> vb) {
    popObj->initialise(vb);
    sort(popObj->population,popObj->population+popObj->pop_size,compare);
    Individual<ZT,FT>v0 = popObj->population[0];
    ZT iter;
    FT prevNorm;
    prevNorm = 0.0;
    while(v0.norm>targetNorm)
    {
        popObj->population = runCrossMut(popObj->dim,k);
        sort(popObj->population,popObj->population+popObj->pop_size,compare);
        cout<<iter<<"\n";
        cout<<v0.norm<<" "<<prevNorm<<"\n";
        v0 = popObj->population[0];
        if(prevNorm!=v0.norm)
        {
             cout<<v0.norm<<" "<<'\n';
             prevNorm = v0.norm;
             iter = 0;
        }
        iter.add_ui(iter,1);
        if(iter==1000)
        {
            cout<<"No. of Iteration exceeds 1000\n";
            cout<<"Restarting Algorithm\n";
            delete []popObj->population;
            return runGA(targetNorm,k,v0);
        }

    }
    return v0.matrix_multiply(v0.x,popObj->get_B(),popObj->dim);
}
// int main()
// {

//     Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>* p1 = new Paper1<Z_NR<mpz_t>,FP_NR<mpfr_t>>("lattices/SVPchallenge",BKZ_DEFAULT,GSO_DEFAULT,10,FT_DEFAULT);
//     Individual<Z_NR<mpz_t>,FP_NR<mpfr_t>>v0;
//     auto start = high_resolution_clock::now();

//     Z_NR<mpz_t>* ans=p1->runGA(2000,10,v0);
//     for(int i=0;i<60;i++)cout<<ans[i]<<" ";
//     cout<<"\n";
//     auto end = high_resolution_clock::now();
//     // Calculate duration and print in milliseconds
//     auto duration = duration_cast<milliseconds>(end - start);
//     cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;
//     cout<<endl;
// }
