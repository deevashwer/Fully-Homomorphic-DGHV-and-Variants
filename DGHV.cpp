#include <fstream>
#include "DGHV_utils.h"

/*
 * For the above (commented) values of parameters, we have a security level of 42 bits
 Currently, I've used rand() (cryptographically insecure), seeded with system time, for this implementation.
 */

using namespace std;

void print(GEN x){
    cout << GENtostr(x) << endl;
    return;
}

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
        GEN ct = generate_secondary_error();
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

    /*
    // Encrypts an integer
    GEN encrypt(GEN m){
        //printf("Checkpoint: Encrypt");
        //cout<<GENtostr(m)<<endl;
        m = decimal_to_binary_rev(m);
        long n = lg(m) - 1;
        GEN ct = cgetg(n + 1, t_VEC);
        for(int i = 0; i < n; i++)
            gel(ct, i + 1) = encrypt_bit(sk, pk, gel(m, i + 1));
        return ct;
    }*/

    // Decrypts a bit
    GEN decrypt_bit(GEN ct){
        GEN pt = Fp_red(ct, sk);
        pt = Fp_red(pt, gen_2);
        return pt;
    }

    /*
    // Decrypts an integer
    GEN decrypt(GEN ct){
        int n = lg(ct) - 1;
        GEN pt = cgetg(n + 1, t_VEC);
        for(int i = 0; i < n; i++)
            gel(pt, i + 1) = decrypt_bit(sk, gel(ct, i + 1));
        return pt;
    }*/

    // Multiply gate for 2 encrypted bits
    GEN multiply_bit(GEN ct_1, GEN ct_2){
        GEN result = gmul(ct_1, ct_2);
        for(int i = 0; i < gamma; i++){
            result = Fp_red(result, pk[gamma + tau - i - 1]);
        }
        //result = Fp_red(result, pk[0]);
        return result;
    }

    // Addition gate for 2 encrypted bits
    GEN add_bit(GEN ct_1, GEN ct_2){
        return gmod(gadd(ct_2, ct_1), pk[0]);
    }

    /*
    // Addition gate for 2 encrypted integers
    GEN addition_gate(GEN ct_1, GEN ct_2, int n){
        if (n == -1)
            n = lg(ct_1) - 1;
        int m = lg(ct_2) - 1;
        //cout << n << m << endl;
        GEN result = cgetg(max(n, m) + 2, t_VEC);
        GEN carry, temp_carry;
        if(n >= m){
            gel(result, 1) = add_bit(gel(ct_1, 1), gel(ct_2, 1), gel(pk, 1));
            carry =  multiply_bit(gel(ct_1, 1), gel(ct_2, 1), pk);
            //cout << GENtostr(decrypt_bit(sk, gel(result, 1))) << " " << GENtostr(decrypt_bit(sk, carry)) << endl;
            for(int i = 1; i < m; i++){
                gel(result, i + 1) = add_bit(gel(ct_1, i + 1), gel(ct_2, i + 1), gel(pk, 1));
                //cout << "R=A+B= " << GENtostr(decrypt_bit(sk, gel(result, i + 1))) << endl;
                temp_carry = multiply_bit(gel(result, i + 1), carry, pk);
                //cout << "C_1=R*C= " << GENtostr(decrypt_bit(sk, temp_carry)) << endl;
                gel(result, i + 1) = add_bit(gel(result, i + 1), carry, gel(pk, 1));
                //cout << "S=R+C= " << GENtostr(decrypt_bit(sk, gel(result, i + 1))) << endl;
                carry =  multiply_bit(gel(ct_1, i + 1), gel(ct_2, i + 1), pk);
                //cout << "C_2=A*B= " << GENtostr(decrypt_bit(sk, carry)) << endl;
                carry = add_bit(carry, temp_carry, gel(pk, 1));
                //cout << "C_1+C_2= " << GENtostr(decrypt_bit(sk, carry)) << endl;
            }
            //cout << "Checkpoint" << endl;
            for(int i = m; i < n; i++){
                gel(result, i + 1) = add_bit(gel(ct_1, i + 1), carry, gel(pk, 1));
                carry = multiply_bit(gel(ct_1, i + 1), carry, pk);
            }
            gel(result, n + 1) = carry;
        }
        else{
            gel(result, 1) = add_bit(gel(ct_1, 1), gel(ct_2, 1), gel(pk, 1));
            carry =  multiply_bit(gel(ct_1, 1), gel(ct_2, 1), pk);
            for(int i = 1; i < n; i++){
                gel(result, i + 1) = add_bit(gel(ct_1, i + 1), gel(ct_2, i + 1), gel(pk, 1));
                temp_carry = multiply_bit(gel(result, i + 1), carry, pk);
                gel(result, i + 1) = add_bit(gel(result, i + 1), carry, gel(pk, 1));
                carry =  multiply_bit(gel(ct_1, i + 1), gel(ct_2, i + 1), pk);
                carry = add_bit(carry, temp_carry, gel(pk, 1));
            }
            for(int i = n; i < m; i++){
                gel(result, i + 1) = add_bit(gel(ct_2, i + 1), carry, gel(pk, 1));
                carry = multiply_bit(gel(ct_2, i + 1), carry, pk);
            }
            gel(result, m + 1) = carry;
        }
        //cout << "Checkpoint" << endl;
        return result;
    }

    // Shifts ct right by i, for conventional multiplication
    GEN multiplication_utility(GEN result, GEN ct, int i){
        int n = lg(ct) - 1;
        GEN temp = cgetg(n + i + 1, t_VEC);
        for(int j = 0; j < i; j++)
            gel(temp, j + 1) = gen_0;
        for(int j = 0; j < n; j++)
            gel(temp, j + 1 + i) = gel(ct, j + 1);
        //cout << "temp=" << GENtostr(decrypt(sk, temp)) << " ";
        result = addition_gate(result, temp, pk, n + i);
        return result;
    }

    // Multiplication gate for two encrypted integers
    GEN multiplication_gate(GEN ct_1, GEN ct_2){
        int n = lg(ct_1) - 1, m = lg(ct_2) - 1;
        //cout << GENtostr(decrypt(sk, ct_1)) << endl;
        //cout << GENtostr(decrypt(sk, ct_2)) << endl;
        GEN result = cgetg(m + n + 1, t_VEC);
        GEN temp = cgetg(n + 1, t_VEC);
        for(int i = 0; i < m + n; i++)
            gel(result, i + 1) = gen_0;
        for(int j = 0; j < n; j++)
            gel(result, j + 1) = multiply_bit(gel(ct_1, j + 1), gel(ct_2, 1), pk);
        //cout << GENtostr(decrypt(sk, result)) << endl;
        for(int i = 1; i < m; i++){
            for(int j = 0; j < n; j++)
                gel(temp, j + 1) = multiply_bit(gel(ct_1, j + 1), gel(ct_2, i + 1), pk);
            result = multiplication_utility(result, temp, pk, sk, i);
            //cout << "result=" << GENtostr(decrypt(sk, result)) << endl;
        }
        return result;
    }*/
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
        result.value = pkc->add_bit(this->value, ct_2.value);
        result.degree = max(this->degree, ct_2.degree);
        return result;
    }

    ciphertext operator*(const ciphertext& ct_2){
        ciphertext result(pkc);
        result.value = pkc->multiply_bit(this->value, ct_2.value);
        result.degree = this->degree + ct_2.degree;
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
    /*ciphertext ct_1(&pkc, stoi(1));
    ciphertext ct_2(&pkc, stoi(0));
    ciphertext ct_3 = ct_1 * ct_2;
    print(ct_3.decrypt());
    cout << ct_3.degree << endl;*/
    ciphertext ct[100];
    ciphertext result(&pkc, stoi(1));
    for(int i = 0; i < 100; i++){
        ct[i].initialize(&pkc, stoi(1));
        result = result * ct[i];
        if(result.decrypt() == gen_0)
            cout << i << endl;
    }
    print(result.decrypt());
    cout << result.degree << endl;
    pari_close();

    return 0;
}