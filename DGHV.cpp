#include <fstream>
#include "DGHV_utils.h"

//#define BITS_IN_LONG 2

/*
 * For the above (commented) values of parameters, we have a security level of 42 bits
 Currently, I've used rand() (cryptographically insecure), seeded with system time, for this implementation.
 */

using namespace std;

//class ciphertext;

class binary_real{
public:    
    GEN decimal;
    int precision;
    vector<bool> value;

    binary_real(){};

    void initialize(GEN num, GEN den, int Precision){
        precision = Precision;
        GEN remainder;
        decimal = gfloor(gdiv(num, den));
        //print(quotient);
        remainder = gmod(num, den);
        //print(remainder);
        vector<bool>::iterator it;

        remainder = gmul(remainder, gen_2);//, remainder);
        //cout << "Checkpoint" << endl;
        for(int i = 0; i < precision; i++){
            while(gcmp(remainder, den) == -1 && i < precision){
                remainder = gmul(gen_2, remainder);
                value.push_back(false);
                i++;
            }
            if(gcmp(remainder, den) == 1){
                value.push_back(true);
                gsubz(remainder, den, remainder);
                remainder = gmul(remainder, gen_2);
            }
            else if(gcmp(remainder, den) == 0){
                value.push_back(true);
                gsubz(remainder, den, remainder);
            }
            else if(gcmp(remainder, gen_0) == 0)
                value.push_back(false);
        }
        
    }

    void print_real(){
        vector<bool>::iterator it;
        cout << GENtostr(decimal) << ".";
        for(it = value.begin(); it != value.end(); it++){
            cout << *it;
        }
        cout << "" << endl;
    }
};

class cryptosystem{
public:
    GEN sk;
    GEN pk[tau + gamma];

    // Constructor
    cryptosystem(){
        pari_sp ltop_super;
        ifstream sk_file("secret_key.txt");
        if(!sk_file) {
            //printf("Error: Secret Key!");
            ltop_super = avma;
            generate_secret_key();
            sk = gerepilecopy(ltop_super, sk);
            sk_file.close();
            ofstream sk_file("secret_key.txt");
            sk_file << GENtostr(sk);
            sk_file.close();
        }
        else{
            string secret_key;
            sk_file >> secret_key;
            sk = gp_read_str(secret_key.c_str());
            //cout << GENtostr(sk) << endl;
            sk_file.close();
        }

        FILE *pk_file = fopen("public_key.txt", "r");
        if(!pk_file){
            //printf("Error: Public Key");
            fclose(pk_file);
            generate_public_key();
            ofstream pk_file("public_key.txt");
            for(int i = 0; i < tau + gamma; i++)
                pk_file << GENtostr(pk[i]) << endl;
            pk_file.close();
        }
        else{
            //pk = cgetg(tau + 1 + gamma, t_VEC);
            for(int i = 0; i < tau + gamma; i++){
                pk[i] = gp_read_stream(pk_file);
            }
            fclose(pk_file);
        }
    }

    // Outputs secret key, an odd integer in the range (2^(eta - 1), 2^eta). 
    GEN generate_secret_key(){
        pari_sp ltop, lbot;
        sk = gshift(stoi(1), eta - 1);
        gettimeofday(&tv, NULL);
        srand(tv.tv_usec + tv.tv_sec*1000000);
        for(int i = eta - 33; i >= 0; i -= 32){
            if(i == eta - 33){
                unsigned temp = rand();
                while((temp/(1 << 31)) > 0)
                    temp = rand();
                sk = gshift(stoi(temp), i);
                continue;
            }
            ltop = avma;
            GEN temp = gshift(stoi(rand()), i);
            lbot = avma;
            gaddz(sk, temp, sk);
            gerepile(ltop, lbot, sk);
        }
        int least_significant_bits;
        do{
            least_significant_bits = rand();
        }while((least_significant_bits & 1) == 0);
        gaddz(sk, stoi(least_significant_bits), sk);
    }

    /*
    Outputs public key
    Public key has tau elements of the form x = q*p + r, q lies in the range (0, 2^(gamma-eta))
    It has gamma elements x_i s.t x_i = q_i*p + r, q_i lies in the range (2^(gamma+i-1-eta), 2^(gamma+i-eta))
    */
    GEN generate_public_key(){
        pari_sp ltop_super, ltop;
        //pk = cgetg(tau + 1 + gamma, t_VEC);
        int index = 0;
        vector<int> even;
        ltop_super = avma;
        GEN max = gen_0;
        ltop = avma;
        do{
            even.clear();
            max = gen_0;
            index = 0;
            for(int i = 0; i < tau; i++) {
                pk[i] = generate_x(sk);
                //cout << i << endl;
                if(mpcmp(max, pk[i]) < 0){
                    if(mpodd(pk[i]) == 0 || mpodd(Fp_red(pk[i], sk)) == 1){
                        max = pk[i];
                        even.push_back(i);
                        //cout << "even " << i << endl;
                        continue;
                    }
                    max = pk[i];
                    index = i;
                    even.clear();
                    //cout << "index " << index << endl;
                }
            }
        }while(even.size() != 0);
        //cout << "Checkpoint" << endl;
        for(int i = 1; i < gamma + 1; i++){
            pk[i + tau - 1] = gmul(generate_x_i(sk, gamma + i), gen_2);
            cout << "Checkpoint " << i << endl;
        }
        GEN temp = pk[0];
        pk[0] = max;
        pk[index] = temp;
        //pk = gerepilecopy(ltop_super, pk);
    }

    // Encrypts a single bit
    GEN encrypt_bit(GEN bit){
        bool included[tau - 1];
        gettimeofday(&tv, NULL);
        //srand(tv.tv_usec + tv.tv_sec*1000000);
        setrand(stoi(tv.tv_usec + tv.tv_sec*1000000));
        for(int i = 0; i < tau; i++){
            if(mpodd(getrand()) == 1)
                included[i] = true;
            else
                included[i] = false;
        }
        GEN ct = generate_random(sigma);
        gmulz(gen_2, ct, ct);
        gaddz(ct, bit, ct);
        //cout << "Checkpoint" << endl;
        for(int i = 1; i < tau; i++){
            if(included[i - 1]) {
                ct = gadd(ct, gmul(gen_2, pk[i]));
            }
            //ct = Fp_red(ct, gel(pk, 1));
        }
        //cout << GENtostr(ct) << endl;
        ct = Fp_red(ct, pk[0]);
        return ct;
    }

    // Decrypts a bit
    GEN decrypt_bit(GEN ct){
        GEN pt = Fp_red(ct, sk);
        pt = Fp_red(pt, gen_2);
        return pt;
    }

    // AND gate for 2 encrypted bits
    GEN AND_GATE(GEN ct_1, GEN ct_2){
        GEN result = gmul(ct_1, ct_2);
        for(int i = 0; i < gamma; i++){
            result = Fp_red(result, pk[gamma + tau - i - 1]);
        }
        //result = Fp_red(result, pk[0]);
        return result;
    }

    // XOR gate for 2 encrypted bits
    GEN XOR_GATE(GEN ct_1, GEN ct_2){
        return gmod(gadd(ct_2, ct_1), pk[0]);
    }

    // OR gate for 2 encrypted bits
    GEN OR_GATE(GEN ct_1,GEN ct_2){
        return XOR_GATE(AND_GATE(XOR_GATE(ct_1, gen_1), XOR_GATE(ct_2, gen_1)), gen_1);
    }

    // NOT gate for an encrypted bit
    GEN NOT_GATE(GEN ct_1){
        return XOR_GATE(ct_1, gen_1);
    }

    void recrypt_util(GEN ct){
        GEN x_p = gshift(gen_1, kappa);
        x_p = diviiround(x_p, sk);
        sparse_matrix S;
        generate_sparse_matrix(S, x_p);
        binary_real z_i[Theta];
        for(int i = 0; i < Theta; i++){
            z_i[i].initialize(gmul(ct, gel(S.Theta_vector, i + 1)), gshift(gen_1, kappa), 7);
            //z_i[i].print_real();
        }

        //return x_p;
    }
};

class ciphertext{
public:
    GEN value;
    int degree;
    cryptosystem* pkc;

    ciphertext(){};

    ciphertext(cryptosystem* PKC){
        pkc = PKC;
    }

    ciphertext(cryptosystem* PKC, GEN m){
        value = PKC->encrypt_bit(m);
        degree = 1;
        pkc = PKC;
    }

    GEN decrypt(){
        return pkc->decrypt_bit(value);
    }

    ciphertext operator+(const ciphertext& ct_2){
        ciphertext result(pkc);
        result.value = pkc->XOR_GATE(this->value, ct_2.value);
        result.degree = max(this->degree, ct_2.degree);
        return result;
    }

    ciphertext operator*(const ciphertext& ct_2){
        ciphertext result(pkc);
        result.value = pkc->AND_GATE(this->value, ct_2.value);
        result.degree = this->degree + ct_2.degree;
        return result;
    }

    ciphertext operator^(const ciphertext& ct_2){
        ciphertext result(pkc);
        result.value = pkc->OR_GATE(this->value, ct_2.value);
        result.degree = this->degree + ct_2.degree;
        return result;
    }

    ciphertext operator~(){
        ciphertext result(pkc);
        result.value = pkc->NOT_GATE(this->value);
        result.degree = this->degree;
        return result;
    }

    void initialize(cryptosystem* PKC, GEN m){
        value = PKC->encrypt_bit(m);
        degree = 1;
        pkc = PKC;
    }
};

int main(){
    pari_init(6000000000, 2);
    cryptosystem pkc;
    ciphertext ct(&pkc, stoi(0));
    //print(ct.value);
    pkc.recrypt_util(ct.value);
    //binary_real x(stoi(1), stoi(4), 10);
    //GEN x = rdivii(gen_1, pkc.sk, 8);
    //setexpo(x, 6);
    //print(x);
    pari_close();

    return 0;
}