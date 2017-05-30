#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <vector>
#include <gmp.h>
#include <gmpxx.h>

#define rho 2//26
#define sigma 5//42
#define eta 3000
#define gamma 9216//147456
#define tau 100//158
#define Theta 20
#define theta 4
#define kappa (gamma + 3)
#define n 6

#define RAND_MAX (1 << rho + 1)

struct timeval tv;
mpz_t zero, one, two;

using namespace std;

void sort_utility(mpz_t array[], int l, int m, int r){
    //cout << l << " " << m << " " << r << endl;
    int i, j, k;
    int n_1 = m - l + 1, n_2 = r - m;
    mpz_t L[n_1], R[n_2];
    for(i = 0; i < n_1; i++){
        mpz_init(L[i]);
        mpz_set(L[i], array[l + i]);
    }
    for (j = 0; j < n_2; j++){
        mpz_init(R[j]);
        mpz_set(R[j], array[m + 1 + j]);
    }
    i = 0, j = 0, k = l;
    while (i < n_1 && j < n_2){
        if (mpz_cmp(L[i], R[j]) < 1)
            mpz_set(array[k++], L[i++]);
        else
            mpz_set(array[k++], R[j++]);
    }
    
    while (i < n_1)
        mpz_set(array[k++], L[i++]);
    
    while (j < n_2)
        mpz_set(array[k++], R[j++]);
    
    for(i = 0; i < n_1; i++)
        mpz_clear(L[i]);
    for (j = 0; j < n_2; j++)
        mpz_clear(R[j]);
}

void sort_huge_numbers(mpz_t array[], int l, int r){
    if(l < r){
        int m = (l + r)/2;
        sort_huge_numbers(array, l, m);
        sort_huge_numbers(array, m + 1, r);
        sort_utility(array, l, m, r);
    }
    return;
}

void generate_random(mpz_t x, int bit_size){
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec + tv.tv_sec*1000000);
    mpz_t tmp;
    mpz_init(tmp);
    if(bit_size < 33){
        mpz_set_ui(tmp, rand());
        int temp = (1 << bit_size);
        mpz_mod_ui(x, tmp, temp);
        mpz_clear(tmp);
        return;
    }
    for(int i = bit_size - 32; i >= 0; i -= 32){
        mpz_set_ui(tmp, rand());
        mpz_mul_2exp(tmp, tmp, i);
        mpz_add(x, x, tmp);
    }
    mpz_add_ui(x, x, rand());
    mpz_clear(tmp);
}

// Outputs integers of the form x = sk*(q) + r, where q is a random integer in the range (0, 2^(gamma-eta)) and r is the primary error
void generate_x(mpz_t x, mpz_t sk){
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec + tv.tv_sec*1000000);
    mpz_t q, r;
    mpz_inits(q, r, NULL);
    generate_random(q, gamma - eta);
    generate_random(r, rho);
    mpz_mul(x, q, sk);
    mpz_add(x, r, x);
    mpz_clears(q, r, NULL);
}

// Outputs x_i s.t x_i = q_i*p + r, q_i lies in the range (2^(gamma+length-1-eta), 2^(gamma+length-eta))
void generate_x_i(mpz_t x_i, mpz_t sk, int length){
    mpz_t tmp, q, r;
    mpz_inits(tmp, q, r, NULL);
    mpz_mul_2exp(q, one, length - 1);
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec + tv.tv_sec*1000000);
    for(int i = length - 32; i >= 0; i -= 32){
        //GEN temp = gshift(stoi(rand()), i);
        mpz_set_ui(tmp, rand());
        mpz_mul_2exp(tmp, tmp, i);
        mpz_add(q, tmp, q);
    }
    mpz_add_ui(q, q, rand());
    generate_random(r, rho);
    mpz_mul(x_i, q, sk);
    mpz_add(x_i, r, x_i);
    mpz_clears(tmp, q, r, NULL);
}

void generate_sparse_matrix(mpz_t u_i[], bool modified_sk[], mpz_t x_p){
    srand(tv.tv_usec + tv.tv_sec*1000000);
    mpz_t theta_vector[theta];
    mpz_t sum;
    mpz_init(sum);
    for(int i = 0; i < theta - 1; i++){
        mpz_init(theta_vector[i]);
        generate_random(theta_vector[i], kappa - eta);
        mpz_add(sum, sum, theta_vector[i]);
    }
    mpz_init(theta_vector[theta - 1]);
    mpz_set(theta_vector[theta - 1], x_p);
    sort_huge_numbers(theta_vector, 0, theta - 1);
    for(int i = theta - 1; i > 0; i--)
        mpz_sub(theta_vector[i], theta_vector[i], theta_vector[i - 1]);
    for(int i = 0; i < Theta; i++)
        modified_sk[i] = false;
    int temp;
    for(int i = 0; i < theta; i++){
        temp = rand() % Theta;
        while(modified_sk[temp] == true)
            temp = rand() % Theta;
        mpz_set(u_i[temp], theta_vector[i]);
        modified_sk[temp] = true;
    }
    for(int i = 0; i < Theta; i++)
        if(modified_sk[i] == false)
            generate_random(u_i[i], kappa - eta);
    for(int i = 0; i < theta; i++)
        mpz_clear(theta_vector[i]);
    return;
}
