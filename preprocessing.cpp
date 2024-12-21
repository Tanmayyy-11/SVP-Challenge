#include "preprocessing.h"


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

template<class ZT,class FT>preProcessing<ZT,FT>::preProcessing(const char *input_filename,int flags_bkz,int flags_gso,int prec,FloatType float_type)
{
    ZZ_mat<mpz_t>A;
    int status = read_file(A, input_filename);
    dim = A.get_rows();
    if(status==0)
    {
        // int status = bkz_reduction(A, 10, flags_bkz, float_type,prec);
        MatGSO<ZT, FT> M(A, U, UT, flags_gso);
        M.update_gso();
        writeMatrixToFile(A,"finaltest_A.txt");
        B = std::make_unique<std::unique_ptr<ZT[]>[]>(dim);
        Bstar = std::make_unique<std::unique_ptr<FT[]>[]>(dim);
        mu = std::make_unique<std::unique_ptr<FT[]>[]>(dim);
        for(int i = 0;i<dim;i++)
        {
            Bstar[i] = std::make_unique<FT[]>(dim); 
            mu[i] = std::make_unique<FT[]>(dim);
            B[i] = std::make_unique<ZT[]>(dim);
        }
        cout<<"Mu Below:"<<endl;
        for(int i = 0;i<dim;i++)
        {
            for(int j = 0;j<dim;j++)
            {
                M.get_gram(Bstar[i][j],i,j);
                M.get_mu(mu[i][j],i,j);
                cout<<mu[i][j]<<" ";
                B[i][j] = A[i][j];
            }
            cout<<endl;
        }
    }
    else
    {
        cout<<"Error in read file\n";
    }
}


template<class ZT,class FT>preProcessing<ZT,FT>::preProcessing(unique_ptr<std::unique_ptr<ZT[]>[]>& old_B,unique_ptr<ZT[]>& vect,int dimension,int flags_gso)
{
    dim = dimension;
    ZZ_mat<mpz_t> A(dim+1, dim); // to copy original basis and new vector
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            A[i][j] = old_B.get()[i][j];  
        }
    }
    for(int i=0;i<dim;i++){
        cout<<vect.get()[i]<<" ";
    }

    cout<<endl;
    
    for(int i=0;i<dim;i++) {
        A[dim][i] = vect.get()[i];
    }
    writeMatrixToFile(A,"matrix_old.txt");
   
    int status = bkz_reduction(A, 2, BKZ_DEFAULT, FT_DEFAULT, 120);
    MatGSO<ZT, FT> M(A, U, UT, flags_gso);
    M.update_gso();
    writeMatrixToFile(A,"updated_A.txt");

    if(status !=0 ){
        cout<<"Error while running LLL"<<endl;
    }

    int zero_row = -1;//store index of zero vector
    int count  =0 ; 
    for (int i=0;i<=dim ;i++){
        bool all_zero = true;
        for(int j=0;j<dim;j++){
            if(A[i][j] != 0){
                all_zero = false;
                break;
            }
        }
        if(all_zero){
            zero_row = i;
            count++;
        }
    }   
    
    
    ZZ_mat<mpz_t> square_A(dim,dim);//get square matrix before applying gso
    for(int i=0;i<=dim;i++){
        for(int j=0;j<dim;j++){
            if(i<zero_row){
                square_A[i][j] = A[i][j];
            }
            else if(i>zero_row){
                square_A[i-1][j] = A[i][j];
            }
        }
    }
    writeMatrixToFile(square_A,"matrix.txt");

    MatGSO<ZT, FT> M2(square_A, U, UT, flags_gso);
    M2.update_gso();
    B = std::make_unique<std::unique_ptr<ZT[]>[]>(dim);
    Bstar = std::make_unique<std::unique_ptr<FT[]>[]>(dim);
    mu = std::make_unique<std::unique_ptr<FT[]>[]>(dim);
    for(int i = 0;i<dim;i++)
    {
        Bstar[i] = std::make_unique<FT[]>(dim); 
        mu[i] = std::make_unique<FT[]>(dim);
        B[i] = std::make_unique<ZT[]>(dim);
    }
    for(int i = 0;i<dim;i++)
    {
        for(int j = 0;j<dim;j++)
        {
            M2.get_gram(Bstar[i][j],i,j);
            M2.get_mu(mu[i][j],i,j);
            B[i][j] = square_A[i][j];
        }
    }

}


template class preProcessing<Z_NR<mpz_t>,FP_NR<mpfr_t>>;
