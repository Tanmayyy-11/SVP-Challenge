#include <iostream>
#include <fplll/fplll.h>
#include <cmath>
#include <fplll.h>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <random>
#include <map>
#include <filesystem>
#include "test_utils.h"

using namespace std;
using namespace fplll;
using namespace std::chrono;

namespace fs = std::filesystem;


void copy_file(const std::string& source_filename, const std::string& destination_filename) {
    // Open the source file for reading
    std::ifstream source_file(source_filename, std::ios::binary);
    if (!source_file.is_open()) {
        std::cerr << "Error: Could not open source file " << source_filename << std::endl;
        return;
    }

    // Open the destination file for writing
    std::ofstream destination_file(destination_filename, std::ios::binary);
    if (!destination_file.is_open()) {
        std::cerr << "Error: Could not open destination file " << destination_filename << std::endl;
        return;
    }

    // Copy the contents from source to destination
    destination_file << source_file.rdbuf();

    // Close both files
    source_file.close();
    destination_file.close();

    std::cout << "File copied from " << source_filename << " to " << destination_filename << std::endl;
}


std::string substringAfterLastSlash(const std::string& s) {
    // Find the position of the last '/'
    size_t pos = s.rfind('/');
    
    // If '/' is found, return the substring after it and remove the last 4 characters
    if (pos != std::string::npos) {
        std::string result = s.substr(pos + 1);
        if (result.length() > 4) {
            return result.substr(0, result.length() - 4);  // Remove the last 4 characters
        }
        return "";  // If the substring is shorter than or equal to 4 characters, return an empty string
    } else {
        return "";  // If there's no '/', return an empty string
    }
}

map <int, mpz_t> m;

void func(string& input_filename2){
    int beta = 40;
    const char* input_filename = input_filename2.c_str();;
    cout<<"processing file "<<input_filename<<endl;
    ZZ_mat<mpz_t>A;
    int status = read_file(A, input_filename);
    int dim = A.get_rows();
    if(dim<65||dim>79) return;

    int flags_bkz = BKZ_DEFAULT;
    FloatType float_type = FT_DEFAULT;
    int prec = 120;
    int flags_gso = GSO_DEFAULT;
    Matrix<Z_NR<mpz_t>>U;
    Matrix<Z_NR<mpz_t>>UT;
    bkz_reduction(A, beta, flags_bkz, float_type,prec);
    MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> M(A, U, UT, flags_gso);
    M.update_gso();

    mpz_t shortest_length;
    mpz_init_set_ui(shortest_length, 0);
    bool is_first = true;
    int shortest_index = -1;

    // Iterate over each row to find the shortest vector
    for (size_t i = 0; i < A.get_rows(); i++) {
        mpz_t length, sum_of_squares;
        mpz_init(length);
        mpz_init_set_ui(sum_of_squares, 0);

        // Compute the Euclidean norm for the i-th row (sum of squares)
        for (size_t j = 0; j < A.get_cols(); j++) {
            mpz_t square;
            mpz_init(square);
            mpz_mul(square, A[i][j].get_data(), A[i][j].get_data()); // A[i][j]^2
            mpz_add(sum_of_squares, sum_of_squares, square);
            mpz_clear(square);
        }

        // Take the square root of the sum to get the norm
        mpz_sqrt(length, sum_of_squares);

        // Check if this is the shortest vector found so far
        if (is_first || mpz_cmp(length, shortest_length) < 0) {
            mpz_set(shortest_length, length);
            shortest_index = i;
            is_first = false;
        }

        mpz_clear(length);
        mpz_clear(sum_of_squares);
    }

    // Output the length of the shortest vector
    gmp_printf("Length of the shortest vector: %Zd\n", shortest_length);
    std::cout << "Index of the shortest vector: " << shortest_index << std::endl;
    std::cout << "dimension of the lattice is " << A.get_rows() << std::endl;

    string newfilename = "/Users/tanmayshrivastav/Desktop/GeneticAlgorithmsForSVP SmartPointers/lattices2_beta50/";
    string lastname = "dim" + to_string(dim);
    lastname += '_';
    // lastname += "beta=";
    // lastname += to_string(beta);
    // lastname += '_';
    lastname += mpz_get_str(NULL, 10, shortest_length);
    lastname += ".txt";
    newfilename += lastname;
    cout<<newfilename<<endl;
    copy_file(input_filename2,newfilename);


    // Clear mpz_t variable
    mpz_clear(shortest_length);

}

void process_all_txt_files(const std::string& folder_path) {
    // Iterate through all files in the directory
    for (const auto& entry : fs::directory_iterator(folder_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".txt") {
            string filename = entry.path().string();
            func(filename); // Call func on each .txt file
        }
    }
}

// int main(){
//     std::string folder_path = "/Users/tanmayshrivastav/Desktop/GeneticAlgorithmsForSVP SmartPointers/lattices"; // Specify the path to your folder
//     process_all_txt_files(folder_path);
//     return 0;

// }