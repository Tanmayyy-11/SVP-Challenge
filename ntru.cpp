// #include <cmath>
// #include <chrono>
// #include <iostream>
// #include<utility>
// #include <cstdlib>
// #include <ctime>
// #include <cstring>
// #include <memory>
// #include<random>
// #include "rng.h"

// using namespace std;

// int mpow(int base, int exp) {

//   int result = 1;
//   while (exp > 0) {
//     if (exp & 1) result = (result * base) ;
//     base = (base * base) ;
//     exp >>= 1;
//   }
//   return result;
// }

// void printvec(vector<int> &v){
//     for(int i=0;i<v.size();i++){
//         cout<<v[i]<<" ";
//     }
//     cout<<endl;
// }

// int get_q_no_error(int d, int p = 3) {
//     /*
//     The function returns the value of q that gives no decryption failure for a variant of NTRU.
//     Input: d - int (order of the group / 3)
//            p - usually 3
//     */
    
//     int value = p * (6 * d + 1);
//     int q = 1 << (static_cast<int>(std::log2(value)) + 1); // Equivalent to 2^(len(bin(value)) - 2) in Python
//     return q;
// }

// vector<int> vector_mul(vector<int> &a, vector<int> &b, int mod){
//     int n = a.size();
//     vector<int> res(n,0);
//     for(int i=0;i<n;i++){
//         for(int j=0;j<n;j++){
//             res[(i+j)%n] = (a[i]*b[j] + res[(i+j)%n] + mod*mod) % mod;
//         }
//     }
//     return res;
// }

// bool is_degree_one(vector<int>& v){
//     int n = v.size();
//     for(int i=1;i<n;i++){
//         if(v[i] != 0) return false;
//     }
//     return v[0] != 0;
// }

// vector<int> multiplyByConstant(vector<int> &b, int constant, int mod){
//     int n = b.size();
//     vector<int> res;
//     for(int i=0;i<n;i++){
//         res.push_back( (b[i]*constant + constant*mod)%mod);
//     }
//     return res;
// }

// int mod_inverse(int a, int m) {
//     int m0 = m;
//     int t, q;
//     int x0 = 0, x1 = 1;

//     // Edge case: if m is 1, the modular inverse doesn't exist
//     if (m == 1) {
//         return 0;
//     }

//     // Apply Extended Euclidean Algorithm
//     while (a > 1) {
//         // q is the quotient
//         q = a / m;
        
//         // t is the remainder, update a and m
//         t = m;
//         m = a % m;
//         a = t;

//         // Update x0 and x1
//         t = x0;
//         x0 = x1 - q * x0;
//         x1 = t;
//     }

//     // Make x1 positive if it's negative
//     if (x1 < 0) {
//         x1 += m0;
//     }

//     return x1;
// }



// void add_vector(vector<int> &a, vector<int>&b, int mod){
//     if(a.size()!= b.size()){
//         cout<<"error in size"<<endl;
//         return;
//     }
//     for(int i=0;i<a.size();i++){
//         a[i] = (a[i]+b[i] + mod)%mod;
//     }
// }

// void sub_vector(vector<int> &a, vector<int>&b, int mod){
//     if(a.size()!= b.size()){
//         cout<<"error in size"<<endl;
//         return;
//     }
//     for(int i=0;i<a.size();i++){
//         a[i] = (a[i]%mod-b[i]%mod + mod)%mod;
//     }
// }

// void multiplyByX(vector<int>& poly) {
//     poly.insert(poly.begin(), 0); // Shift all coefficients by 1 to multiply by X...lose the coeff for x^n if exists
//     poly.pop_back();
// }

// void divideByX(vector<int>& poly) {//delete first element that was 0...maintain length by appending a 0
//     if (!poly.empty()) {
//         poly.erase(poly.begin());
//         poly.push_back(0);
//     }
// }

// int deg(vector<int> &v){
//     for(int i= v.size()-1;i>-1;i--){
//         if(v[i]!=0){
//             return i;
//         }
//     }
// }

// vector<int> roundMultiplication(vector<int> &a, int power, int mod){
    
//     int n = a.size() -1 ;
//     vector<int> result (n+1,0);
//     for(int i=0;i<n+1;i++){
//         result[(i+power)%n] = (a[i] + result[(i+power)%n])% mod;
//     }
//     result[0] = (result[0]+result[n])%mod;
//     result[n] = 0;
//     return result;

// }

// vector<int> minus_roundMultiplication(vector<int> &a, int power, int mod){
    
//     vector<int> result = roundMultiplication(a,power,mod);
//     for(int i=0;i<result.size();i++){
//         result[i] = (-result[i]%mod + mod)%mod;
//     }
//     return result;

// }

// bool is_one(vector<int> &f){//check f(x) = 1 
//     int n = f.size();
//     for(int i=1;i<n;i++){
//         if(f[i]!=0)return false;
//     }
//     return (f[0] == 1);
// }

// bool is_minusone(vector<int> &f, int mod){//check f(x) = 1 
//     int n = f.size();
//     for(int i=1;i<n;i++){
//         if(f[i]!=0)return false;
//     }
//     return (f[0]%mod == (-1 + mod)%mod);
// }


// class Poly{
//     public:
//     int n;               // order
//     vector<int> coeff; // Coefficients
//     Poly(int order) : n(order), coeff(order, 0) {}

//     void sample(vector<bool> &b, int d1, int d2){// d1 number of 1 and d2 number of -1
//         for(int i=0;i<n;i++){
//             for(int j =0 ;j<30;j++){
//                 coeff[i] += (1<<(2+j)) * b[30*i + j];
//             }
//             if(i<d1)coeff[i]++;
//             else if(i<d1+d2)coeff[i]+=2;
//         }

//         sort(coeff.begin(),coeff.end());

//         for(int i=0;i<n;i++){
//             coeff[i]&=3; //mod4
//             if(coeff[i]==2) coeff[i] = -1;
//         }
//     }

//     void multiply(Poly* a, Poly*b, int q){
//         for(int i=0;i<n;i++) coeff[i] =0;
//         for(int i=0;i<n;i++){
//             for(int j=0;j<n;j++){
//                 coeff[(i+j)%n] = (coeff[(i+j)%n] + a->coeff[i]*b->coeff[j])%q;
//             }
//         }
//     }

//     vector<int> inverse2(){
//         vector<int> a (n+1, 0); //initialize a vector of size n+1 is taking quotient x^n - 1 and copy the contents from coeff of 1 to x^n-1
//         for(int i=0;i<n;i++){
//             a[i] = coeff[i];
//         }
//         int k=0;
//         vector<int> b (n+1, 0);//intialize b as 1... all other coeff as 0
//         b[0] = 1;
//         vector<int> c (n+1, 0);//intialize c as 0
//         vector<int>& f = a ;//reference the same vector a
//         vector<int> g (n+1, 0);// represent x*N -1 
//         g[0] = 1;
//         g[n] = 1;

//         bool t = true;

//         while(t){

//             while(f[0] == 0){   //step3
//                 divideByX(f);   //step4
//                 multiplyByX(c);
//                 k++;
//             }

//             if(is_one(f)){
//                 vector<int> ans =  roundMultiplication(b,n-k,2);//step 5
//                 ans.pop_back();
//                 return ans;
//             }

//             if(deg(f) < deg(g)){// step 6
//                 swap(f,g);  //step 7
//                 swap(b,c);
//             }

//             add_vector(f,g,2);
//             add_vector(b,c,2);

//         }

//     }


//     vector<int> inverse3(){
//         int mod = 3;
//         vector<int> ans;
//         vector<int> a (n+1, 0); //initialize a vector of size n+1 is taking quotient x^n - 1 and copy the contents from coeff of 1 to x^n-1
//         for(int i=0;i<n;i++){
//             a[i] = coeff[i];
//         }
//         int k=0;
//         vector<int> b (n+1, 0);//intialize b as 1... all other coeff as 0
//         b[0] = 1;
//         vector<int> c (n+1, 0);//intialize c as 0
//         vector<int>f (a);//copy
//         vector<int> g (n+1, 0);// represent x*N -1 
//         g[0] = 2;
//         g[n] = 1;

//         bool t = true;

//         while(t){

//             while(f[0] == 0){   //step3
//                 divideByX(f);   //step4
//                 multiplyByX(c);
//                 k++;
//             }

//             if(is_one(f)){
//                 vector<int> ans= roundMultiplication(b,n-k,mod);//step 5
//                 ans.pop_back();
//                 return ans;
//             }
//             else if(is_minusone(f,mod)){
//                 vector<int> ans =  minus_roundMultiplication(b,n-k,mod);
//                 ans.pop_back();
//                 return ans;
//             }

//             if(deg(f) < deg(g)){// step 6
//                 swap(f,g);  //step 7
//                 swap(b,c);
//             }

//             if(f[0] == g[0]){
//                 sub_vector(f,g,mod);
//                 sub_vector(b,c,mod);
//             }
//             else{
//                 add_vector(f,g,3);
//                 add_vector(b,c,3);
//             }
//         }
//     }

//     vector<int> inversep(int mod){
//         vector<int> a (n+1, 0); //initialize a vector of size n+1 is taking quotient x^n - 1 and copy the contents from coeff of 1 to x^n-1
//         for(int i=0;i<n;i++){
//             a[i] = coeff[i];
//         }
//         int k=0;
//         vector<int> b (n+1, 0);//intialize b as 1... all other coeff as 0
//         b[0] = 1;
//         vector<int> c (n+1, 0);//intialize c as 0

//         int upperbound = n + deg(a);
//         int counter = 0;

//         vector<int>& f =a;

//         vector<int> g (n+1, 0);// represent x*N -1 
//         g[0] = 2;
//         g[n] = 1;

//         bool t = true;

//         while(t){

//             while(f[0] == 0){   //step3
//                 divideByX(f);   //step4
//                 counter++;
//                 if(counter == upperbound){
//                     cout<<"inverse doesnt exist"<<endl;
//                     return {};
//                 }
//                 multiplyByX(c);
//                 k++;
//             }

//             if(is_degree_one(f)){

//                 int f0_inv = mod_inverse(f[0],mod);
//                 b = multiplyByConstant(b,f0_inv,mod);
//                 vector<int> ans = roundMultiplication(b,n-k,mod);//step 5
//                 ans.pop_back();
//                 return ans;
//             }

//             if(deg(f) < deg(g)){// step 6
//                 swap(f,g);  //step 7
//                 swap(b,c);
//             }

//             int g0_inv = mod_inverse(g[0],mod);
//             int u = f[0]*g0_inv % mod;
//             vector<int> gTemp = multiplyByConstant(g,u,mod);
//             sub_vector(f,gTemp,mod);
//             vector<int> cTemp = multiplyByConstant(c,u,mod);
//             sub_vector(b,cTemp,mod);

//         }

//     }

//     vector<int> inverse_p_power_r(int p, int r){

//         vector<int> a (n, 0); //initialize a vector of size n+1 is taking quotient x^n - 1 and copy the contents from coeff of 1 to x^n-1
//         for(int i=0;i<n;i++){
//             a[i] = coeff[i];
//         }

//         vector<int> b = this->inverse2();
//         printvec(b);
//         if(!b.size()) return b;
//         int q = p;
//         int bound  = mpow(p,r);
//         cout<<"here";
//         while(q<bound){
//             q = q*q;
//             vector<int> t = vector_mul(a,b,q);
//             t[0] = 2 - t[0]; //check
//             for(int i=1;i<t.size();i++) t[i] = -t[i];
//             b = vector_mul(b,t,q);
//         }
//         return b;
//     }
// };

// class NTRU{
//     public:
//     Poly *f;
//     Poly *g;
//     Poly *h;
//     int p,q,n,d;
//     unsigned long long sample_size;
//     vector<vector<int>> basis;
//     unsigned char *seed;
//     vector<bool> b;


//     NTRU(int order, int modulusP, int modulusQ, int dVal): p(modulusP), q(modulusQ), n(order), d(dVal),basis(2 * order, std::vector<int>(2 * order, 0)) 
//     {
//         f = new Poly(order);
//         g = new Poly(order);
//         h = new Poly(order);
//         sample_size = ceil((30 * n) / 8);
//         seed = new unsigned char[sample_size];
//         randombytes(seed,sample_size);
//         for (size_t i = 0; i < sample_size; ++i) {
//             unsigned char byte = seed[i];
//             for (int bit = 7; bit >= 0; --bit) {
//                 b.push_back((byte >> bit) & 1);
//             }
//         }
//     }

//     void sample(){
//         f->sample(b,d+1,d);
//         g->sample(b,d,d);
//     }

// };

// // int main(){
// //     int n = 11;
// //     int p =3;
// //     int d = n/3;
// //     int q = 64;
// //     int x = 6;
    
// //     // NTRU* n1 = new NTRU(11,p,q,d);
// //     // vector<int> f_inverse_p;
// //     // vector<int> f_inverse_q;
// //     // while(true){
// //     //     n1->sample();
// //     //     f_inverse_p = n1->f->inversep(p);
// //     //     f_inverse_q = n1->f->inverse_p_power_r(2,1);
// //     //     vector<int> modp_ans = vector_mul(n1->f->coeff, f_inverse_p, p);
// //     //     vector<int> modq_ans = vector_mul(n1->f->coeff, f_inverse_q, q);
// //     //     if(f_inverse_p.size()>0 && f_inverse_q.size()>0 && is_one(modp_ans) && is_one(modq_ans)){
// //     //         break;
// //     //     }
// //     // }
// //     // printvec(n1->f->coeff);
// //     // printvec(f_inverse_p);
// //     // printvec(f_inverse_q);
// //     // printvec(n1->g->coeff);


// //     // vector<int> mul = vector_mul(n1->f->coeff, f_inverse_p, p);
// //     // printvec(mul);

// //     // for(int i=0;i<n1->f->coeff.size();i++){
// //     //     cout<<n1->f->coeff[i]<<" ";
// //     // }
// //     // cout<<endl;
// //     // for(int i=0;i<n1->g->coeff.size();i++){
// //     //     cout<<n1->g->coeff[i]<<" ";
// //     // }

// //     // Poly* a =new Poly(3);
// //     // a->coeff = {0,1,0};
// //     // vector<int> ans = a->inverse2();
// //     // for(int i=0;i<ans.size();i++){
// //     //     cout<<ans[i]<<" ";
// //     // }
// //     // cout<<endl;

// //     Poly* b =new Poly(11);
// //     b->coeff = {1, -1, 0, 0, 1, -1, 1, 1, -1, 0, 0};
// //     vector<int> ans2 = b->inverse_p_power_r(2,6);
// //     // vector<int> ans2 = b->inverse2();
// //     for(int i=0;i<ans2.size();i++){
// //         cout<<ans2[i]<<" ";
// //     }
// //     cout<<endl;


  

// //     return 0;
// // }


