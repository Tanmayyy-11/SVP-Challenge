#include "Individual.h"
template<class ZT,class FT>
unique_ptr<ZT[]> Individual<ZT,FT>::matrix_multiply(ZT* vec, unique_ptr<ZT[]>* mat, int dim) {
    unique_ptr<ZT[]> result = make_unique<ZT[]>(dim);
    for (int i = 0; i < dim; ++i) {
        result[i] = 0;
        for (int j = 0; j < dim; ++j) {
            ZT temp; 
            temp.mul(vec[j], mat[j][i]);
            result[i].add(result[i], temp);
        }
    }
    return result;
}

template<class ZT,class FT>
unique_ptr<ZT[]> Individual<ZT,FT>::YtoX(ZT* y,unique_ptr<FT[]>* mu, int dim) {
    unique_ptr<ZT[]> x = make_unique<ZT[]>(dim);
    // cout<<"ti below:"<<endl;
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
        // cout<<i<<" "<<t_z<<endl;
        x[i].sub(y[i],t_z);
    }
    return x;
}
template<class ZT,class FT>
unique_ptr<ZT[]> Individual<ZT,FT>::XtoY(ZT * x,unique_ptr<FT[]>* mu, int dim) {
    unique_ptr<ZT[]>y = make_unique<ZT[]>(dim);
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
        y[i].add(x[i],t_z);
    }
    return y;
}
template<class ZT,class FT>
Individual<ZT,FT>::Individual(int dim, unique_ptr<FT[]>* mu, FT* alpha, unique_ptr<ZT[]>* B, unique_ptr<FT[]>* Bstar) {

    this->x = make_unique<ZT[]>(dim);
    this->y = make_unique<ZT[]>(dim);
    this->dim = dim;
    unique_ptr<ZT[]> alpha_r=make_unique<ZT[]>(dim);
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
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dist(0, 1);//to flip sign of bit
        if (dist(gen) == 1) {
            mpz_neg(y[i].get_data(), y[i].get_data()); 
        }
        cout<<y[i]<<" ";
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
    unique_ptr<ZT[]> vect = (matrix_multiply(x.get(), B, dim));
    this->norm = get_norm(vect.get(), dim);
    // delete[] vect;
}
template<class ZT,class FT>
Individual<ZT,FT>::~Individual() 
{
}
template<class ZT,class FT>
Individual<ZT,FT>::Individual()
{
    // cout<<"Hello\n";
    this->x = nullptr;
    this->y = nullptr;
    this->norm = 0.0;
}
template<class ZT,class FT>
Individual<ZT,FT>::Individual(int dim)
{
    this->x = make_unique<ZT[]>(dim);
    this->y = make_unique<ZT[]>(dim);
    this->norm = 0.0;
}
// template<class ZT,class FT>
// Individual<ZT,FT>::Individual(const Individual &t)
// {
//     x = make_unique<ZT[]>(dim);
//     y = make_unique<ZT[]>(dim);
//     norm = t.norm;
//     for(int i = 0;i<dim;i++)
//     {
//         x[i] = t.x[i];
//         y[i] = t.y[i];

//     }
    
// }
template<class ZT, class FT>
Individual<ZT, FT>::Individual(const Individual &t) {
    // Copy the dimension and norm from the original object
    // cout<<"1\n";
    this->dim = t.dim;
    this->norm = t.norm;

    // Deep copy the unique_ptr<ZT[]> x
    if (t.x) {
        // Allocate new memory for x in the copied object
        // cout<<dim<<"\n";
        this->x = std::make_unique<ZT[]>(dim);
        // Copy the elements from the original object's x array
        // std::copy(t.x.get(), t.x.get() + dim, this->x.get());
        // cout<<"X copied\n";
        // cout<<t.x[0]<<"sd\n";
        
        for(int i = 0;i<dim;i++)
        {
            // cout<<t.x[i]<<" ";
            x[i] = t.x[i];
        }
        // cout<<"X copied\n";

        
    } else {
        // If the original x was null, keep the new one null as well
        this->x = nullptr;
    }

    // Deep copy the unique_ptr<ZT[]> y
    if (t.y) {
        // Allocate new memory for y in the copied object
        // cout<<dim<<"\n";
        this->y = std::make_unique<ZT[]>(dim);
        // Copy the elements from the original object's y array
        // std::copy(t.y.get(), t.y.get() + dim, this->y.get());
        for(int i = 0;i<dim;i++)y[i] = t.y[i];

    } else {
        // If the original y was null, keep the new one null as well
        this->y = nullptr;
    }
    // cout<<"Copy Done\n";
}
template<class ZT, class FT>
Individual<ZT, FT>::Individual(Individual &&t) noexcept : dim(t.dim), norm(t.norm), x(std::move(t.x)), y(std::move(t.y)) {
    // cout<<"2\n";
    // Leave t in a valid state
    t.dim = 0;
    t.norm = FT();  // Reset the norm
}
template<class ZT, class FT>
Individual<ZT, FT>& Individual<ZT, FT>::operator=(const Individual &t) {
    // cout<<"3\n";
    if (this != &t) {
        dim = t.dim;
        norm = t.norm;
        if (t.x) {
            x = std::make_unique<ZT[]>(dim);
            for(int i = 0;i<dim;i++)x[i] = t.x[i];
        }
        if (t.y) {
            y = std::make_unique<ZT[]>(dim);
            for(int i = 0;i<dim;i++)y[i] = t.y[i];
        }
    }
    // cout<<"done\n";
    return *this;
}
template<class ZT, class FT>
Individual<ZT, FT>& Individual<ZT, FT>::operator=(Individual &&t) noexcept {
    // cout<<"4\n";
    if (this != &t) {
        dim = t.dim;
        norm = t.norm;
        x = std::move(t.x);
        y = std::move(t.y);
        
        // Leave t in a valid state
        t.dim = 0;
        t.norm = FT();
    }
    return *this;
}
template class Individual<Z_NR<mpz_t>,FP_NR<mpfr_t>>;
