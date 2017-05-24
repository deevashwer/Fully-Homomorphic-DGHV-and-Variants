#include <fstream>
#include "DGHV_utils.h"

//#define BITS_IN_LONG 2

/*
 * For the above (commented) values of parameters, we have a security level of 42 bits
 Currently, I've used rand() (cryptographically insecure), seeded with system time, for this implementation.
 */

using namespace std;

struct encrypted_sparse_matrix{
    GEN encrypted_secret_key[Theta];
    GEN encrypted_z[Theta + 1][8];
};

class binary_real{
public:    
    GEN decimal;
    int precision;
    vector<bool> value;

    binary_real(){};

    binary_real(GEN num, GEN den, int Precision){
        precision = Precision;
        GEN remainder;
        decimal = gfloor(gdiv(num, den));
        //print(quotient);
        remainder = gmod(num, den);
        //print(remainder);
        vector<bool>::iterator it;

        remainder = gmul(remainder, gen_2);//, remainder);
        //cout << "Checkpoint" << endl;
        int i = 0;
        while(i < precision){
            if(gcmp(remainder, den) == -1){
                remainder = gmul(gen_2, remainder);
                value.push_back(false);
            }
            else if(gcmp(remainder, den) == 1){
                value.push_back(true);
                gsubz(remainder, den, remainder);
                remainder = gmul(remainder, gen_2);
            }
            else if(gcmp(remainder, den) == 0){
                value.push_back(true);
                gsubz(remainder, den, remainder);
            }
            else if(gcmp(remainder, gen_0) == 0){
                value.push_back(false);
            }
            i++;
        }
        return;
    }

    void initialize(GEN num, GEN den, int Precision){
        precision = Precision;
        //print(num);
        //print(den);
        GEN remainder;
        decimal = gfloor(gdiv(num, den));
        //print(decimal);
        remainder = gmod(num, den);
        //print(remainder);
        vector<bool>::iterator it;

        remainder = gmul(remainder, gen_2);//, remainder);
        //cout << "Checkpoint" << endl;
        int i = 0;
        while(i < precision){
            if(gcmp(remainder, den) == -1){
                remainder = gmul(gen_2, remainder);
                value.push_back(false);
            }
            else if(gcmp(remainder, den) == 1){
                value.push_back(true);
                gsubz(remainder, den, remainder);
                remainder = gmul(remainder, gen_2);
            }
            else if(gcmp(remainder, den) == 0){
                value.push_back(true);
                gsubz(remainder, den, remainder);
            }
            else if(gcmp(remainder, gen_0) == 0){
                value.push_back(false);
            }
            i++;
        }
        return;
    }

    void custom_setup(GEN dec, int prec, vector<bool> val){
        decimal = dec;
        precision = prec;
        value = val;
    }

    binary_real operator+(const binary_real& x){
        binary_real result;
        GEN dec = gadd(this->decimal, x.decimal);
        int prec = this->precision;
        vector<bool> bits, val;
        bool carry = false;
        
        /*
        Assuming all the truncated bits to be 1. This gets a initial true carry for the additon of the least significant bit.
        This adds an error of size 1/16(Theta) to the sum.
        
        bool carry = true;*/
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
            gaddz(dec, gen_1, dec);
        for(iterator = bits.rbegin(); iterator != bits.rend(); iterator++)
            val.push_back(*iterator);
        result.custom_setup(dec, prec, val);
        return result;
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
        srand(tv.tv_usec + tv.tv_sec*1000000);
        for(int i = 0; i < tau; i++)
            if((rand() & 1) == 1)
                included[i] = true;
            else
                included[i] = false;

        GEN ct = generate_random(sigma);
        gmulz(gen_2, ct, ct);
        gaddz(ct, bit, ct);

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

    GEN squashed_decryption(GEN ct){
        GEN x_p = gshift(gen_1, kappa);
        x_p = diviiround(x_p, sk);
        sparse_matrix S;
        generate_sparse_matrix(S, x_p);
        int n = ceil(log2(theta)) + 3;
        binary_real z_i[Theta];
        for(int i = 0; i < Theta; i++){
            z_i[i].initialize(gmul(ct, gel(S.Theta_vector, i + 1)), gshift(gen_1, kappa), n);
        }
        binary_real result;
        result.initialize(gen_0, gen_1, ceil(log2(theta)) + 3);
        for(int i = 0; i < Theta; i++){
            if(gel(S.modified_secret_key, i + 1) == gen_1)
                result = result + z_i[i];
        }
        //result.print_real();
        if(result.value.at(0) == true)
            gaddz(result.decimal, gen_1, result.decimal);
        return Fp_red(gsub(ct, result.decimal), gen_2);
    }

    encrypted_sparse_matrix recrypt_util(GEN ct, cryptosystem* PKC){
        GEN x_p = gshift(gen_1, kappa);
        x_p = diviiround(x_p, sk);
        sparse_matrix S;
        generate_sparse_matrix(S, x_p);
        //for(int i = 0; i < Theta; i++)
        //    cout << GENtostr(gel(S.Theta_vector, i + 1)) << endl;
        int n = ceil(log2(theta)) + 3;
        binary_real z_i[Theta];
        for(int i = 0; i < Theta; i++){
            z_i[i].initialize(gmul(ct, gel(S.Theta_vector, i + 1)), gshift(gen_1, kappa), n);
            //z_i[i].print_real();
        }
        encrypted_sparse_matrix a_i;
        //a_i.encrypted_secret_key = cgetg(Theta + 1, t_VEC);
        //encrypted_secret_key = new ciphertext[Theta];
        for(int i = 0; i < Theta; i++)
            a_i.encrypted_secret_key[i] = PKC->encrypt_bit(gel(S.modified_secret_key, i + 1));
        //a_i.encrypted_z = cgetg(Theta + 1, t_MAT);
        //for(int i = 0; i < Theta; i++)
        //    gel(a_i.encrypted_z, i + 1) = cgetg(n + 1, t_COL);
        for(int i = 0; i < Theta; i++){
            for(int j = 0; j < n + 1; j++){
                if(j == n){
                    if(mpodd(z_i[i].decimal) == 1){
                        a_i.encrypted_z[i][j] = PKC->encrypt_bit(gen_1);
                    }
                    else{
                        a_i.encrypted_z[i][j] = PKC->encrypt_bit(gen_0);
                    }
                }
                else{
                    //cout << z_i[i].value.at(n - j - 1);
                    a_i.encrypted_z[i][j] = PKC->encrypt_bit(stoi(z_i[i].value.at(n - j - 1)));
                }
            }
        }
        /*for(int i = 0; i < Theta; i++){
            for(int j = 0; j < n + 1; j++){
                cout << GENtostr(PKC->decrypt_bit(a_i.encrypted_z[i][j])) << " ";
            }
            cout << "" << endl;
        }*/
        return a_i;         
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

    void custom_setup(GEN val, int deg, cryptosystem* PKC){
        value = val;
        degree = deg;
        pkc = PKC;
    }

    void recrypt(cryptosystem* PKC){
        int n = ceil(log2(theta)) + 3;
        ciphertext encrypted_secret_key[Theta];
        ciphertext encrypted_z[Theta][n + 1];
        pari_sp ltop;
        //ltop = avma;
        encrypted_sparse_matrix a_i = pkc->recrypt_util(value, PKC);
        //encrypted_sparse_matrix = gerepilecopy(ltop, encrypted_sparse_matrix);
        ciphertext encrypted_a_i[Theta][n + 1];
        for(int i = 0; i < Theta; i++)
            encrypted_secret_key[i].custom_setup(a_i.encrypted_secret_key[i], 1, PKC);
        cout << "Checkpoint" << endl;
        for(int i = 0; i < Theta; i++){
            for(int j = 0; j < n + 1; j++){
                encrypted_z[i][j].custom_setup(a_i.encrypted_z[i][j], 1, PKC);
                encrypted_a_i[i][j] = encrypted_secret_key[i] * encrypted_z[i][j];
                cout << i << " " << j << endl;
            }
        }
        cout << "Checkpoint" << endl;
        ciphertext encrypted_ciphertext_bit;
        if(mpodd(value) == 1)
            encrypted_ciphertext_bit.initialize(PKC, gen_1);
        else
            encrypted_ciphertext_bit.initialize(PKC, gen_0);
        ciphertext dp[n + 2][n + 2];
        dp[0][0].initialize(PKC, gen_1);
        cout << "Checkpoint" << endl;
        for(int i = 1; i < n + 2; i++){
            dp[i][0].initialize(PKC, gen_0);
            dp[0][i].initialize(PKC, gen_1);
        }
        cout << "Checkpoint" << endl;
        for(int j = 1; j < n + 2; j++)
            for(int i = 1; i < n + 2; i++){
                dp[i][j] = (encrypted_a_i[i - 1][j - 1] * dp[i - 1][j - 1]) + (dp[i][j - 1]);
                cout << i << ", " << j << " " << dp[i][j].degree << endl;
            }
        cout << "Checkpoint" << endl;
        cout << dp[n + 1][n + 1].degree << endl;

    }
};

int main(){
    pari_init(60000000000, 2);
    cryptosystem pkc;
    ciphertext ct(&pkc, stoi(1));
    //print(pkc.squashed_decryption(ct.value));
    //binary_real b(ct.value, pkc.sk, 3000);
    //b.print_real();
    ct.recrypt(&pkc);
    //binary_real b(stoi(3), stoi(7), 6);
    //b.print_real();
    pari_close();

    return 0;
}