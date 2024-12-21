#include "repres.h"
using namespace std;
using namespace fplll;
void writeMatrixTo(ZZ_mat<mpz_t>& matrix, const std::string& filename) {
    std::ofstream outFile(filename);  // Open file for writing

    if (!outFile) {  // Check if the file opened successfully
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    int dim = matrix.get_rows();  // Get the number of rows
    outFile << dim << " " << dim << std::endl;  // Write dimensions

    for (int i=0;i<dim;i++) {  // Iterate through each row
        for (int j=0;j<dim-1;j++) {  // Iterate through each element in the row
            outFile << matrix[i][j] << " ";  // Write the element
        }
        outFile << std::endl;  // New line after each row
    }

    outFile.close();  // Close the file
    std::cout << "Matrix written to " << filename << std::endl;
}


template<class ZT,class FT>repres<ZT,FT>::repres(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type)
{

    preprocess = make_unique<preProcessing<ZT,FT>>(input_filename, flags_bkz, flags_gso, prec, float_type);
    dim = preprocess->dim;
    pop_size = 2*dim;
    FT normB0 = get_norm((preprocess->Bstar[0]).get(),dim);
    FT log2 = log(2);
    alpha=make_unique<FT[]>(dim);
    length=make_unique<ZT[]>(dim);
    totLength = 0;
    for(int i = 0;i<dim;i++)
    {
        alpha[i] = normB0/get_norm((preprocess->Bstar[i]).get(),dim);
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
unique_ptr<ZT[]>  repres<ZT,FT>::decode(bool *chromosome) {
    unique_ptr<ZT[]> y = make_unique<ZT[]>(dim);
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
unique_ptr<bool[]>  repres<ZT,FT>::encode(ZT* y, ZT totalLength) {
    unique_ptr<bool[]>  chromosome = make_unique<bool[]>((int)totalLength.get_d());
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

    population = make_unique<Individual<ZT,FT>[]>(pop_size);
    if(v0.x==NULL)
    {
        population[0].y = make_unique<ZT[]>(dim) ;
        population[0].y[0] = 1;
        for(int i = 1;i<dim;i++)population[0].y[i] = 0;

        population[0].x = unique_ptr<ZT[]>(population[0].YtoX(population[0].y.get(),preprocess->mu.get(),dim));
        unique_ptr<ZT[]> vect = unique_ptr<ZT[]>(population[0].matrix_multiply(population[0].x.get(), (preprocess->B).get(),dim));
        population[0].norm = population[0].get_norm(vect.get(),dim);
        cout<<population[0].norm<<" first vector norm"<<endl;
        population[0].dim = dim;
    }
    else{ 
        population[0] = v0;

        //update the basis by applying LLL
        //either get Nelite to the next generation but don't forget to update y sparse representation
        //or randomly generate new population with respect to the new basis
        for(int i=0;i<dim;i++){
            cout<<population[0].x.get()[i]<<" ";
        }
        cout<<"below xs"<<endl;
         ZZ_mat<mpz_t>B_mat(dim,dim);
         for(int i=0;i<dim;i++){
            for(int j=0;j<dim;j++){
                B_mat[i][j] = (preprocess->B).get()[i][j];
            }
         }
         writeMatrixTo(B_mat, "Output_B.txt");
        unique_ptr<ZT[]> vect = unique_ptr<ZT[]>(population[0].matrix_multiply(population[0].x.get(), (preprocess->B).get(),dim));
        unique_ptr<preProcessing<ZT,FT>> new_preprocess = make_unique<preProcessing<ZT,FT>>(preprocess->B,vect,dim,GSO_DEFAULT);
        // preprocess = make_unique<preProcessing<ZT,FT>>(preprocess->B,vect,dim,GSO_DEFAULT);
        preprocess = move(new_preprocess);
        
        FT normB0 = get_norm((preprocess->Bstar[0]).get(),dim);
        cout<<normB0<<endl;
        FT log2 = log(2);
        totLength = 0;
        for(int i = 0;i<dim;i++)
        {
            alpha[i] = normB0/get_norm((preprocess->Bstar[i]).get(),dim);
            alpha[i]=sqrt(alpha[i]);
            length[i] .set_f(floor(log(alpha[i])/log2)+FT(2));
            totLength.add(totLength,length[i]);
        }
        cout<<"here"<<endl;





    //     for(int i = 0;i<dim;i++)
    //     {
    //         population[0].x[i] = v0.x[i];
    //         population[0].y[i] = v0.y[i];
    //     }
    //    population[0].norm =  v0.norm ;
    //    population[0].dim = v0.dim;
    }
    
    for (int i = 1; i < pop_size; i++) {
        
     population[i] =  Individual<ZT, FT>(dim, (preprocess->mu).get(), alpha.get(), (preprocess->B).get(), (preprocess->Bstar).get());
     if(get_norm(population[i].x.get(),dim) == 0.0){
        cout<<" i = "<<i<<endl;
        i--;
     }
    }
}
template<class ZT,class FT>unique_ptr<ZT[]>* repres<ZT,FT>::get_B()
{
    return preprocess->B.get();
}
template<class ZT,class FT> unique_ptr<FT[]>* repres<ZT,FT>::get_mu()
{
    return preprocess->mu.get();
}

template class repres<Z_NR<mpz_t>,FP_NR<mpfr_t>>;
