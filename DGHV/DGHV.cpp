#include <gmp.h>
//#include <stdio.h>
#include <fstream>
#include <cmath>
#include "DGHV_utils.h"

//int n = ceil(log2(theta)) + 3;

using namespace std;

class binary_real{
public:
    int decimal;
    int precision;
    vector<bool> value;
    
    binary_real(){};
    
    binary_real(mpz_t num, mpz_t den, int Precision){
        precision = Precision;
        mpz_t remainder, quotient;
        mpz_inits(remainder, quotient, NULL);
        mpz_fdiv_q(quotient, num, den);
        if(mpz_odd_p(quotient) == 1)
            decimal = 1;
        else
            decimal = 0;
        mpz_fdiv_r(remainder, num, den);
        vector<bool>::iterator it;
        
        mpz_mul(remainder, remainder, two);
        int i = 0;
        while(i < precision){
            if(mpz_cmp(remainder, den) == -1){
                mpz_mul(remainder, remainder, two);
                value.push_back(false);
            }
            else if(mpz_cmp(remainder, den) == 1){
                value.push_back(true);
                mpz_sub(remainder, remainder, den);
                mpz_mul(remainder, remainder, two);
            }
            else if(mpz_cmp(remainder, den) == 0){
                value.push_back(true);
                mpz_sub(remainder, remainder, den);
            }
            else if(mpz_cmp(remainder, zero) == 0){
                value.push_back(false);
            }
            i++;
        }
        mpz_clears(remainder, quotient, NULL);
        return;
    }
    
    void initialize(mpz_t num, mpz_t den, int Precision){
        precision = Precision;
        mpz_t remainder, quotient;
        mpz_inits(remainder, quotient, NULL);
        mpz_fdiv_q(quotient, num, den);
        if(mpz_odd_p(quotient) == 1)
            decimal = 1;
        else
            decimal = 0;
        mpz_fdiv_r(remainder, num, den);
        vector<bool>::iterator it;
        
        mpz_mul(remainder, remainder, two);
        int i = 0;
        while(i < precision){
            if(mpz_cmp(remainder, den) == -1){
                mpz_mul(remainder, remainder, two);
                value.push_back(false);
            }
            else if(mpz_cmp(remainder, den) == 1){
                value.push_back(true);
                mpz_sub(remainder, remainder, den);
                mpz_mul(remainder, remainder, two);
            }
            else if(mpz_cmp(remainder, den) == 0){
                value.push_back(true);
                mpz_sub(remainder, remainder, den);
            }
            else if(mpz_cmp(remainder, zero) == 0){
                value.push_back(false);
            }
            i++;
        }
        mpz_clears(remainder, quotient, NULL);
        return;
    }
    
    void custom_setup(int dec, int prec, vector<bool> val){
        decimal = dec;
        precision = prec;
        value = val;
    }
    
    binary_real operator+(const binary_real& x){
        binary_real result;
        int dec = (this->decimal + x.decimal) % 2;
        int prec = this->precision;
        vector<bool> bits, val;
        //bool carry = false;
        
        /*
         Assuming all the truncated bits to be 1. This gets a initial true carry for the additon of the least significant bit.
         This adds an error of size 1/16(Theta) to the sum.
         */
         bool carry = true;
        vector<bool>::reverse_iterator iterator;
        for(int i = 0; i < prec; i++){
            if(this->value.at(precision - i - 1) + x.value.at(precision - i - 1) + carry == 3){
                bits.push_back(true);
                carry = true;
            }
            else if(this->value.at(precision - i - 1) + x.value.at(precision - i - 1) + carry == 2){
                bits.push_back(false);
                carry = true;
            }
            else{
                bits.push_back(this->value.at(precision - i - 1) + x.value.at(precision - i - 1) + carry);
                carry = false;
            }
        }
        if(carry == true)
            dec = (dec + 1) % 2;
        for(iterator = bits.rbegin(); iterator != bits.rend(); iterator++)
            val.push_back(*iterator);
        result.custom_setup(dec, prec, val);
        return result;
    }
    
    void print_real(){
        vector<bool>::iterator it;
        cout << decimal << ".";
        for(it = value.begin(); it != value.end(); it++){
            cout << *it;
        }
        cout << "" << endl;
    }
};

class cryptosystem{
public:
    mpz_t sk;
    mpz_t pk[tau + gamma];
    
    cryptosystem(){
        mpz_init(sk);
        for(int i = 0; i < tau + gamma; i++)
            mpz_init(pk[i]);
        ifstream sk_file("secret_key.txt");
        if(!sk_file) {
            generate_secret_key();
            sk_file.close();
            ofstream sk_file("secret_key.txt");
            sk_file << mpz_get_str(NULL, 10, sk);
            sk_file.close();
        }
        else{
            string secret_key;
            sk_file >> secret_key;
            mpz_set_str(sk, secret_key.c_str(), 10);
            sk_file.close();
        }
        cout << "Secret Key Generated." << endl;
        
        ifstream pk_file("public_key.txt");
        if(!pk_file){
            pk_file.close();
            generate_public_key();
            ofstream pk_file("public_key.txt");
            for(int i = 0; i < tau + gamma; i++)
                pk_file << mpz_get_str(NULL, 10, pk[i]) << endl;
            pk_file.close();
        }
        else{
            string public_key;
            for(int i = 0; i < tau + gamma; i++){
                pk_file >> public_key;
                mpz_set_str(pk[i], public_key.c_str(), 10);
                //cout << i << " " << pk[i] << endl;
            }
            pk_file.close();
        }
        cout << "Public Key Generated." << endl;
    }
    
    // Outputs secret key, an odd integer in the range (2^(eta - 1), 2^eta).
    void generate_secret_key(){
        mpz_t tmp;
        mpz_init(tmp);
        mpz_mul_2exp(sk, one, eta - 1);
        gettimeofday(&tv, NULL);
        srand(tv.tv_usec + tv.tv_sec*1000000);
        for(int i = eta - 32; i >= 0; i -= 32){
            mpz_set_ui(tmp, rand());
            mpz_mul_2exp(tmp, tmp, i);
            mpz_add(sk, tmp, sk);
        }
        mpz_add_ui(sk, sk, rand());
        mpz_mul(sk, sk, two);
        mpz_add(sk, sk, one);
        mpz_clear(tmp);
    }
    
    /*
     Outputs public key
     Public key has tau elements of the form x = q*p + r, q lies in the range (0, 2^(gamma-eta))
     It has gamma elements x_i s.t x_i = q_i*p + r, q_i lies in the range (2^(gamma+i-1-eta), 2^(gamma+i-eta))
     */
    void generate_public_key(){
        int index = 0;
        vector<int> even;
        mpz_t max, tmp;
        mpz_inits(max, tmp, NULL);
        do{
            even.clear();
            mpz_set_ui(max, 0);
            index = 0;
            for(int i = 0; i < tau; i++) {
                mpz_init(pk[i]);
                generate_x(pk[i], sk);
                //cout << i << endl;
                if(mpz_cmp(max, pk[i]) < 0){
                    mpz_mod(tmp, pk[i], sk);
                    if(mpz_odd_p(pk[i]) == 0 || mpz_odd_p(tmp) == 1){
                        mpz_set(max, pk[i]);
                        even.push_back(i);
                        //cout << "even " << i << endl;
                        continue;
                    }
                    mpz_set(max, pk[i]);
                    index = i;
                    even.clear();
                    //cout << "index " << index << endl;
                }
            }
        }while(even.size() != 0);
        //cout << "Checkpoint " << endl;
        for(int i = 1; i < gamma + 1; i++){
            generate_x_i(pk[tau + i - 1], sk, gamma - eta + i - 1);
            mpz_mul(pk[tau + i - 1], pk[tau + i - 1], two);
            //cout << "Checkpoint " << i << endl;
        }
        mpz_clears(tmp, max, NULL);
        swap(pk[0], pk[index]);
        /*for(int i = 1; i < tau; i++)
            if(mpz_cmp(pk[0], pk[i]) < 0)
                cout << i << " Dayum!" << endl;*/
    }
    
    // Encrypts a single bit
    void encrypt_bit(mpz_t ct, mpz_t bit){
        bool included[tau - 1];
        gettimeofday(&tv, NULL);
        srand(tv.tv_usec + tv.tv_sec*1000000);
        for(int i = 0; i < tau; i++)
            if((rand() & 1) == 1)
                included[i] = true;
            else
                included[i] = false;
        mpz_t tmp;
        mpz_init(tmp);
        generate_random(ct, sigma);
        mpz_mul(ct, ct, two);
        mpz_add(ct, ct, bit);
        for(int i = 1; i < tau; i++)
            if(included[i - 1]) {
                mpz_mul(tmp, pk[i], two);
                mpz_add(ct, ct, tmp);
                //mpz_mod(tmp, ct, sk);
                //cout << i << " " << tmp << endl;
            }
        mpz_mod(ct, ct, pk[0]);
        mpz_clear(tmp);
    }
    
    // Decrypts a bit
    void decrypt_bit(mpz_t m, mpz_t ct){
        mpz_mod(m, ct, sk);
        mpz_mod(m, m, two);
    }
    
    // AND gate for 2 encrypted bits
    void AND_GATE(mpz_t result, mpz_t ct_1, mpz_t ct_2){
        mpz_mul(result, ct_1, ct_2);
        for(int i = 0; i < gamma; i++)
            mpz_mod(result, result, pk[gamma + tau - i - 1]);
        mpz_mod(result, result, pk[0]);
    }
    
    // XOR gate for 2 encrypted bits
    void XOR_GATE(mpz_t result, mpz_t ct_1, mpz_t ct_2){
        mpz_add(result, ct_1, ct_2);
        mpz_mod(result, result, pk[0]);
    }
    
    // OR gate for 2 encrypted bits
    void OR_GATE(mpz_t result, mpz_t ct_1, mpz_t ct_2){
        mpz_t tmp;
        mpz_init(tmp);
        XOR_GATE(result, ct_1, one);
        XOR_GATE(tmp, one, ct_2);
        AND_GATE(result, result, tmp);
        XOR_GATE(result, result, one);
        mpz_clear(tmp);
    }
    
    // NOT gate for an encrypted bit
    void NOT_GATE(mpz_t result, mpz_t ct_1){
        XOR_GATE(result, ct_1, one);
    }
    
    void recrypt_util(mpz_t encrypted_sk[], mpz_t encrypted_z[][n + 1], mpz_t ct, cryptosystem* PKC){
        mpz_t x_p;
        mpz_init(x_p);
        mpz_mul_2exp(x_p, one, kappa);
        mpz_fdiv_q(x_p, x_p, sk);
        mpz_t u_i[Theta];
        for(int i = 0; i < Theta; i++)
            mpz_init(u_i[i]);
        bool modified_sk[Theta];
        generate_sparse_matrix(u_i, modified_sk, x_p);
        binary_real z_i[Theta];
        mpz_t num, den;
        mpz_inits(num, den, NULL);
        mpz_mul_2exp(den, one, kappa);
        for(int i = 0; i < Theta; i++){
            mpz_mul(num, u_i[i], ct);
            z_i[i].initialize(num, den, n);
            /*if(modified_sk[i] == true){
            //    cout << i << " ";
                z_i[i].print_real();
            }*/
        }
        mpz_clears(num, den, NULL);
        for(int i = 0; i < Theta; i++)
            if(modified_sk[i] == true)
                PKC->encrypt_bit(encrypted_sk[i], one);
            else
                PKC->encrypt_bit(encrypted_sk[i], zero);
        for(int i = 0; i < Theta; i++)
            for(int j = 0; j < n + 1; j++)
                if(j == 0)
                    if((z_i[i].decimal) == 1)
                        PKC->encrypt_bit(encrypted_z[i][j], one);
                    else
                        PKC->encrypt_bit(encrypted_z[i][j], zero);
                else
                    if((z_i[i].value.at(j - 1)) == true)
                        PKC->encrypt_bit(encrypted_z[i][j], one);
                    else
                        PKC->encrypt_bit(encrypted_z[i][j], zero);
    }
};

void two_for_three_trick(mpz_t p, mpz_t q, mpz_t a, mpz_t b, mpz_t c, cryptosystem* pkc){
    mpz_t temp_1, temp_2, temp_3;
    mpz_inits(temp_1, temp_2, temp_3, NULL);
    pkc->XOR_GATE(temp_2, a, b);
    pkc->XOR_GATE(temp_2, temp_2, c);
    pkc->AND_GATE(temp_1, a, b);
    pkc->XOR_GATE(temp_3, a, b);
    pkc->AND_GATE(temp_3, temp_3, c);
    pkc->XOR_GATE(temp_1, temp_1, temp_3);
    mpz_set(p, temp_1);
    mpz_set(q, temp_2);
    return;
}

class ciphertext{
public:
    mpz_t value;
    int degree;
    cryptosystem* pkc;
    
    ciphertext(){
        mpz_init(value);
    };
    
    ciphertext(cryptosystem* PKC){
        mpz_init(value);
        pkc = PKC;
    }
    
    ciphertext(cryptosystem* PKC, mpz_t m){
        mpz_init(value);
        PKC->encrypt_bit(value, m);
        degree = 1;
        pkc = PKC;
    }
    
    void decrypt(mpz_t m){
        pkc->decrypt_bit(m, value);
    }
    
    void initialize(cryptosystem* PKC, mpz_t m){
        PKC->encrypt_bit(value, m);
        degree = 1;
        pkc = PKC;
    }
    
    void custom_setup(mpz_t val, int deg, cryptosystem* PKC){
        mpz_init(value);
        mpz_set(value, val);
        degree = deg;
        pkc = PKC;
        //mpz_clear(val);
    }
    
    void print(){
        mpz_t m;
        mpz_init(m);
        decrypt(m);
        cout << "Degree: " << degree << ", Decrypted Bit: " << m << endl;
        mpz_mod(m, value, pkc->sk);
        cout << "Error: " << m << endl;
        mpz_clear(m);
    }
    
    ciphertext operator+(ciphertext& ct_2){
        ciphertext result(pkc);
        pkc->XOR_GATE(result.value, this->value, ct_2.value);
        result.degree = max(this->degree, ct_2.degree);
        return result;
    }
    
    ciphertext operator*(ciphertext& ct_2){
        ciphertext result(pkc);
        pkc->AND_GATE(result.value, this->value, ct_2.value);
        result.degree = this->degree + ct_2.degree;
        return result;
    }
    
    ciphertext operator^(ciphertext& ct_2){
        ciphertext result(pkc);
        pkc->OR_GATE(result.value, this->value, ct_2.value);
        result.degree = this->degree + ct_2.degree;
        return result;
    }
    
    ciphertext operator~(){
        ciphertext result(pkc);
        pkc->NOT_GATE(result.value, this->value);
        result.degree = this->degree;
        return result;
    }
    
    void clean(){
        mpz_clear(value);
    }
    
    void recrypt(cryptosystem* PKC){
        mpz_t encrypted_sk[Theta];
        mpz_t encrypted_z[Theta][n + 1];
        for(int i = 0; i < Theta; i++){
            for(int j = 0; j < n + 1; j++)
                mpz_init(encrypted_z[i][j]);
            mpz_init(encrypted_sk[i]);
        }
        pkc->recrypt_util(encrypted_sk, encrypted_z, value, PKC);
        mpz_t temp;
        mpz_init(temp);
        ciphertext a[Theta][n + 1];
        for(int i = 0; i < Theta; i++){
            //cout << i << " ";
            for(int j = 0; j < n + 1; j++){
                PKC->AND_GATE(encrypted_z[i][j], encrypted_z[i][j], encrypted_sk[i]);
                a[i][j].custom_setup(encrypted_z[i][j], 2, PKC);
                //PKC->decrypt_bit(temp, a[i][j].value);
                //cout << temp << " ";
                //mpz_clear(encrypted_z[i][j]);
                //a[i][j].print();
            }
            //cout << endl;
            mpz_clear(encrypted_sk[i]);
        }
        ciphertext dp[(int) pow(2, n - 4) + 1][Theta + 1];
        ciphertext W[n + 1][n + 1];
        //int dummy = 1;
        for(int i = 1; i < n - 1; i++)
            dp[i][0].custom_setup(zero, 1, PKC);
        for(int i = 0; i < Theta + 1; i++)
            dp[0][i].custom_setup(one, 1, PKC);
        for(int k = 0; k < n + 1; k++){
            for(int i = 1; i < pow(2, n - 4) + 1; i++){
                for(int j = 1; j < Theta + 1; j++){
                    //cout << k << " " << i << " " << j << endl;
                    dp[i][j] = a[j - 1][k] * dp[i - 1][j - 1];// + dp[i][j - 1];
                    dp[i][j] = dp[i][j] + dp[i][j - 1];
                    /*if(j == Theta && i == dummy){
                        cout << k << " " << i << " " << j << " ";
                        dp[i][j].print();
                        dummy *= 2;
                    }*/
                }
            }
            //dummy = 1;
            if(k < n - 3){
                for(int i = 0; i < k + 1; i++)
                    W[k][i].custom_setup(dp[(int) pow(2, k - i)][Theta].value, dp[(int) pow(2, k - i)][Theta].degree, PKC);
            }
            else{
                for(int i = 0; i < n - 3; i++)
                    W[k][k - n + 4 + i].custom_setup(dp[(int) pow(2, n - 4 - i)][Theta].value, dp[(int) pow(2, n - 4 - i)][Theta].degree, PKC);
            }
            //for(int i = 0; i < Theta; i++)
                //a[i][k].clean();
        }
        /*for(int i = 0; i < n + 1; i++){
            for(int j = 0; j < n + 1; j++){
                PKC->decrypt_bit(temp, W[i][j].value);
                cout << temp << " ";
            }
            cout << endl;
        }*/
        int k = 0, l = 0, size = n + 1;
        while(size > 2){
            k = 0, l = 0;
            while(size > (k + 2)){
                PKC->XOR_GATE(W[l + 1][0].value, W[k][0].value, W[k + 1][0].value);
                PKC->XOR_GATE(W[l + 1][0].value, W[l + 1][0].value, W[k + 2][0].value);
                //W[l + 1][0] = W[k][0] + W[k + 1][0] + W[k + 2][0];
                for(int j = 1; j < n + 1; j++)
                    two_for_three_trick(W[l][j - 1].value, W[l + 1][j].value, W[k][j].value, W[k + 1][j].value, W[k + 2][j].value, PKC);
                mpz_set(W[l][n].value, zero);
                l += 2;
                k += 3;
            }
            if((k + 2) == size){
                for(int j = 0; j < n + 1; j++){
                    mpz_set(W[l][j].value, W[k][j].value);
                    mpz_set(W[l + 1][j].value, W[k + 1][j].value);
                }
                l += 2;
            }
            else if((k + 1) == size){
                for(int j = 0; j < n + 1; j++)
                    mpz_set(W[l][j].value, W[k][j].value);
                l += 1;
            }
            size = l;
        }
        /*cout << endl;
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < n + 1; j++){
                PKC->decrypt_bit(temp, W[i][j].value);
                cout << temp << " ";
            }
            cout << endl;
        }*/
        
        mpz_t c_p_bit;
        mpz_init(c_p_bit);
        PKC->AND_GATE(c_p_bit, W[0][1].value, W[1][1].value);
        PKC->XOR_GATE(c_p_bit, c_p_bit, W[1][0].value);
        PKC->XOR_GATE(c_p_bit, c_p_bit, W[0][0].value);
        if(mpz_odd_p(value) == 0)
            mpz_add(c_p_bit, c_p_bit, one);
        mpz_set(value, c_p_bit);
        pkc = PKC;
        return;
    }
};

int main(){
    mpz_inits(zero, one, two, NULL);
    mpz_set_ui(zero, 0);
    mpz_set_ui(one, 1);
    mpz_set_ui(two, 2);
    cryptosystem pkc;
    //ciphertext ct(&pkc, one);
    //ct.recrypt(&pkc);
    ciphertext ct[350];
    ciphertext result(&pkc, zero);
    for(int i = 0; i < 350; i++){
        ct[i].initialize(&pkc, one);
        result = (result * ct[i]);
        //result.print();
    }
    result.print();
    result.recrypt(&pkc);
    result.print();
    return 0;
}
