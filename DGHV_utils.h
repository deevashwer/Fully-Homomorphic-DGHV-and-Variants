#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <vector>
#include <gmp.h>
#include <gmpxx.h>

#define lambda 42
#define alpha 42
#define rho 16//26
#define beta 12
#define sigma 42//42
#define eta 1088
#define gamma 160000
#define tau 158//158
#define Theta 144
#define theta 15
#define n 4
#define kappa (/*gamma*/161569 + 2 + n)
#define e 0

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

int bit_size(mpz_t x){
    int temp = 0;
    mpz_t tmp;
    mpz_init(tmp);
    if(mpz_cmp(x, zero) < 0)
        mpz_mul_si(tmp, x, -1);
    else
        mpz_set(tmp, x);
    while(mpz_cmp(tmp, zero) > 0){
        mpz_fdiv_q(tmp, tmp, two);
        temp++;
    }
    mpz_clear(tmp);
    return temp;
}

void mpz_mod_modified(mpz_t rop, mpz_t op1, mpz_t op2){
    mpz_mod(rop, op1, op2);
    mpz_t temp;
    mpz_init(temp);
    mpz_fdiv_q(temp, op2, two);
    if(mpz_cmp(rop, temp) > 0)
        mpz_sub(rop, rop, op2);
    mpz_clear(temp);
    return;
}

void generate_random(mpz_t x, int bit_size, bool include_negative_range, bool seeded, bool full_range){
    gettimeofday(&tv, NULL);
    if(!seeded)
        srand(tv.tv_usec + tv.tv_sec*1000000);
    mpz_set(x, zero);
    if(full_range)
        mpz_mul_2exp(x, one, bit_size - 1);
    mpz_t tmp;
    mpz_init(tmp);
    if(bit_size < 33){
        mpz_set_ui(tmp, rand());
        int temp = (1 << bit_size);
        mpz_mod_ui(x, tmp, temp);
        if(include_negative_range){
            mpz_mul_2exp(tmp, one, bit_size - 1);
            mpz_sub(x, x, tmp);
        }
        mpz_clear(tmp);
        return;
    }
    for(int i = bit_size - 32; i >= 0; i -= 32){
        mpz_set_ui(tmp, rand());
        mpz_mul_2exp(tmp, tmp, i);
        mpz_add(x, x, tmp);
    }
    mpz_add_ui(x, x, rand());
    if(include_negative_range){
        mpz_mul_2exp(tmp, one, bit_size - 1);
        mpz_sub(x, x, tmp);
    }
    mpz_clear(tmp);
    return;
}

// Outputs integers of the form x = sk*(q) + r, where q is a random integer in the range (0, 2^(gamma-eta)) and r is the primary error
void generate_x(mpz_t x, mpz_t sk){
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec + tv.tv_sec*1000000);
    mpz_t q, r;
    mpz_inits(q, r, NULL);
    generate_random(q, gamma - eta, false, false, false);
    generate_random(r, rho, true, false, false);
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
    generate_random(r, rho, true, false, false);
    mpz_mul(x_i, q, sk);
    mpz_add(x_i, r, x_i);
    mpz_clears(tmp, q, r, NULL);
}

int generate_sparse_matrix(mpz_t u_1, bool modified_secret_key[], mpz_t x_p){
    int seed = time(NULL);
    srand(seed);
    mpz_t Theta_vector[Theta];
    //mpz_init(Theta_vector[0]);
    for(int i = 1; i < Theta; i++){
        mpz_init(Theta_vector[i]);
        generate_random(Theta_vector[i], kappa + 1, false, true, false);
    }
    modified_secret_key[0] = true;
    for(int i = 1; i < Theta; i++)
        modified_secret_key[i] = false;
    int count = theta - 1, index;
    srand(tv.tv_usec + tv.tv_sec*1000000);
    while(count > 0){
        index = rand() % Theta;
        if(!modified_secret_key[index]){
            modified_secret_key[index] = true;
            count--;
        }
    }
    mpz_t sum, temp;
    mpz_inits(sum, temp, NULL);
    mpz_set(sum, zero);
    mpz_mul_2exp(temp, one, kappa + 1);
    for(int i = 1; i < Theta; i++)
        if(modified_secret_key[i] == true)
            mpz_add(sum, sum, Theta_vector[i]);
    mpz_mod(sum, sum, temp);
    mpz_sub(u_1, x_p, sum);
    if(mpz_cmp(u_1, zero) < 0)
        mpz_add(u_1, temp, u_1);
    /*mpz_set(sum, u_1);
    for(int i = 1; i < Theta; i++)
        if(modified_secret_key[i] == true)
            mpz_add(sum, sum, Theta_vector[i]);
    mpz_mod(sum, sum, temp);
    cout << mpz_cmp(x_p, sum) << endl;*/
    for(int i = 1; i < Theta; i++)
        mpz_clear(Theta_vector[i]);
    mpz_clears(temp, sum, NULL);
    return seed;
}
