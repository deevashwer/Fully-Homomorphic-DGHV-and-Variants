#include <iostream>
#include <cstdlib>
#include <pari/pari.h>
#include <ctime>
#include <sys/time.h>
#include <vector>

#define rho 10//26
#define sigma 20//42
#define eta 988
#define gamma 9216//147456
#define tau 100//158
#define Theta 150
#define theta 15
#define kappa (gamma + 2)

#define RAND_MAX (1 << rho + 1)

struct timeval tv;

using namespace std;

// Generate errors between (0, 2^rho)
int generate_error(){
    return (rand() % RAND_MAX);// - (RAND_MAX / 2);
}

// Returns a reverse bit vector of integer m
GEN decimal_to_binary_rev(GEN m){
    if(gcmp(m, gen_0) == 0){
        GEN binary_m = cgetg(2, t_VEC);
        gel(binary_m, 1) = gen_0;
        return binary_m;
    }
    vector<bool> bits;
    GEN temp;
    while(mpcmp(m, gen_1) > 0){
        if(mpodd(m) == 1)
            bits.push_back(true);
        else
            bits.push_back(false);
        m = gdivexact(m, gen_2);
    }
    bits.push_back(true);
    GEN binary_m = cgetg(bits.size() + 1, t_VEC);
    int i = 1;
    vector<bool>::iterator it;
    for(it = bits.begin(); it != bits.end(); it++)
        if(*it == true)
            gel(binary_m, i++) = gen_1;
        else
            gel(binary_m, i++) = gen_0;
    bits.clear();
    return binary_m;
}

// Returns secondary error of sigma bits, used at the time of encryption. 
GEN generate_secondary_error(){
    pari_sp ltop, lbot;
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec + tv.tv_sec*1000000);
    GEN r = stoi(rand());
    for(int i = sigma - 33; i >= 0; i -= 32){
        ltop = avma;
        GEN temp = gshift(stoi(rand()), i);
        lbot = avma;
        gaddz(r, temp, r);
        gerepile(ltop, lbot, r);
    }
    return r;
}


// Outputs integers of the form x = sk*(q) + r, where q is a random integer in the range (0, 2^(gamma-eta)) and r is the primary error
GEN generate_x(GEN sk){
    pari_sp ltop, lbot, ltop_super;
    ltop_super =  avma;
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec + tv.tv_sec*1000000);
    GEN q;
    for(int i = gamma - eta; i >= 0; i -= 32){
        if(i == gamma - eta){
            q = gshift(stoi(rand()), i);
            continue;
        }
        ltop = avma;
        GEN temp = gshift(stoi(rand()), i);
        lbot = avma;
        gaddz(q, temp, q);
        gerepile(ltop, lbot, q);
    }
    gaddz(q, stoi(rand()), q);
    GEN r = stoi(generate_error());
    GEN x = gmul(sk, q);
    gaddz(x, r, x);
    //cout << GENtostr(r) << endl;
    //cout << GENtostr(Fp_red(x, sk)) << endl;
    x = gerepilecopy(ltop_super, x);
    return x;
}

// Outputs x_i s.t x_i = q_i*p + r, q_i lies in the range (2^(gamma+length-1-eta), 2^(gamma+length-eta))
GEN generate_x_i(GEN sk, int length){
    pari_sp ltop, lbot, ltop_super;
    ltop_super =  avma;
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec + tv.tv_sec*1000000);
    GEN q;
    for(int i = length - eta; i >= 0; i -= 32){
        if(i == length - eta){
            unsigned temp = rand();
            while((temp/(1 << 31)) > 0)
                temp = rand();
            q = gshift(stoi(temp), i);
            continue;
        }
        ltop = avma;
        GEN temp = gshift(stoi(rand()), i);
        lbot = avma;
        gaddz(q, temp, q);
        gerepile(ltop, lbot, q);
    }
    gaddz(q, stoi(rand()), q);
    GEN r = stoi(generate_error());
    GEN x = gmul(sk, q);
    gaddz(x, r, x);
    //cout << GENtostr(r) << endl;
    //cout << GENtostr(Fp_red(x, sk)) << endl;
    x = gerepilecopy(ltop_super, x);
    return x;
}
