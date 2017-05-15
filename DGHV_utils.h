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
#define kappa (gamma + 3)

#define RAND_MAX (1 << rho + 1)

struct timeval tv;

struct sparse_matrix{
    GEN Theta_vector;
    GEN modified_secret_key;
};

using namespace std;

void print(GEN x){
    cout << GENtostr(x) << endl;
    return;
}

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

// Returns random integer of bit_length bits.
GEN generate_random(int bit_length){
    gettimeofday(&tv, NULL);
    setrand(stoi(tv.tv_usec + tv.tv_sec*1000000));
    GEN r = randomi(gshift(gen_1, bit_length));
    return r;
}


// Outputs integers of the form x = sk*(q) + r, where q is a random integer in the range (0, 2^(gamma-eta)) and r is the primary error
GEN generate_x(GEN sk){
    gettimeofday(&tv, NULL);
    setrand(stoi(tv.tv_usec + tv.tv_sec*1000000));
    GEN q = generate_random(gamma - eta);
    GEN r = generate_random(rho);
    GEN x = gmul(sk, q);
    gaddz(x, r, x);
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

int calculate_bit_length(GEN x){
    int bit_length = 1;
    while(mpcmp(x, gen_1) > 0){
        x = gdivexact(x, gen_2);
        bit_length++;
    }
    return bit_length;
}

void generate_sparse_matrix(sparse_matrix& S, GEN x_p){
    //sparse_matrix S;
    GEN theta_vector = cgetg(theta + 1, t_VEC);
    setrand(stoi(tv.tv_usec + tv.tv_sec*1000000));
    for(int i = 0; i < theta - 1; i++)
        gel(theta_vector, i + 1) = randomi(x_p);
    gel(theta_vector, theta) = x_p;
    theta_vector = sort(theta_vector);
    for(int i = theta; i > 1; i--)
        gsubz(gel(theta_vector, i), gel(theta_vector, i - 1), gel(theta_vector, i));
    S.modified_secret_key = cgetg(Theta + 1, t_VEC);
    for(int i = 0; i < Theta; i++)
        gel(S.modified_secret_key, i + 1) = gen_0;
    S.Theta_vector = cgetg(Theta + 1, t_VEC);
    ulong temp;
    for(int i = 0; i < theta; i++){
        temp = random_Fl(Theta);
        while(mpcmp(gel(S.modified_secret_key, temp + 1), gen_0) != 0)
            temp = random_Fl(Theta);
        gel(S.Theta_vector, temp + 1) = gel(theta_vector, i + 1);
        gel(S.modified_secret_key, temp + 1) = gen_1;
    }
    for(int i = 0; i < Theta; i++)
        if(mpcmp(gel(S.modified_secret_key, i + 1), gen_0) == 0)
            gel(S.Theta_vector, i + 1) = randomi(x_p);
    return;// S;
}
