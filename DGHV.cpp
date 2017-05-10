#include <iostream>
#include <cstdlib>
#include <pari/pari.h>
#include <ctime>
#include <sys/time.h>
#include <vector>
#include <fstream>

#define rho 13//26
#define sigma 26//42
#define eta 250//988
#define gamma 9216//147456
#define tau 100//158
#define Theta 150
#define theta 15
#define kappa (gamma + 2)

#define RAND_MAX (1 << rho + 1)

struct timeval tv;
//GEN total_error = gen_0;

/*
 * For the above mentioned values of parameters, we have a security level of 42 bits
 */

using namespace std;

int generate_error(){
    return (rand() % RAND_MAX);// - (RAND_MAX / 2);
}

GEN decimal_to_binary_rev(GEN m){
    //cout << "Checkpoint " << GENtostr(m) << " " << gcmp(m, stoi(6)) << endl;
    if(gcmp(m, gen_0) == 0){
        GEN binary_m = cgetg(2, t_VEC);
        gel(binary_m, 1) = gen_0;
        return binary_m;
    }
    vector<bool> bits;
    //int count = 5;
    GEN temp;
    while(mpcmp(m, gen_1) > 0){
    //while(gcmp(m, gen_0) > 0 && count-- > 0){
        if(mpodd(m) == 1)
            bits.push_back(true);
        else
            bits.push_back(false);
        m = gdivexact(m, gen_2);
        //cout << "Checkpoint " << GENtostr(m) << " " << gcmp(gadd(m, gen_1), gen_1) << endl;
    }
    //cout << "Checkpoint " << mpcmp(gen_0, stoi(1)) << endl;
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
        //cout << i << "\t" << GENtostr(temp) << endl;
    }
    return r;
}

GEN generate_secret_key(){
    pari_sp ltop, lbot;
    GEN sk = gshift(stoi(1), eta - 1);
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec + tv.tv_sec*1000000);
    for(int i = eta - 33; i >= 0; i -= 32){
        ltop = avma;
        GEN temp = gshift(stoi(rand()), i);
        lbot = avma;
        gaddz(sk, temp, sk);
        gerepile(ltop, lbot, sk);
        //cout << i << "\t" << GENtostr(temp) << endl;
    }
    int least_significant_bits;

    do{
        least_significant_bits = rand();
    }while((least_significant_bits & 1) == 0);
    gaddz(sk, stoi(least_significant_bits), sk);
    return sk;
}

GEN generate_x(GEN sk, int length){
    pari_sp ltop, lbot, ltop_super;
    ltop_super =  avma;
    gettimeofday(&tv, NULL);
    srand(tv.tv_usec + tv.tv_sec*1000000);
    GEN q;
    for(int i = length - eta; i >= 0; i -= 32){
        if(i == length - eta){
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

GEN key_gen(GEN sk){
    pari_sp ltop_super, ltop;
    GEN pk = cgetg(tau + 1 + gamma, t_VEC);
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
            gel(pk, i + 1) = generate_x(sk, gamma);
            //cout << i << endl;
            if(mpcmp(max, gel(pk, i + 1)) < 0){
                if(mpodd(gel(pk, i + 1)) == 0 || mpodd(Fp_red(gel(pk, i + 1), sk)) == 1){
                    max = gel(pk, i + 1);
                    even.push_back(i);
                    //cout << "even " << i << endl;
                    continue;
                }
                max = gel(pk, i + 1);
                index = i;
                even.clear();
                //cout << "index " << index << endl;
            }
        }
    }while(even.size() != 0);
    cout << "Checkpoint" << endl;
    for(int i = 1; i < gamma + 1; i++){
        gel(pk, i + tau) = gmul(generate_x(sk, gamma + i), gen_2);
        cout << "Checkpoint " << i << endl;
    }
    GEN temp = gel(pk, 1);
    gel(pk, 1) = max;
    gel(pk, index + 1) = temp;
    pk = gerepilecopy(ltop_super, pk);
    return pk;
}

GEN encrypt_bit(GEN sk, GEN pk, GEN bit){
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

    for(int i = 1; i < tau; i++){
        if(included[i - 1]) {
            ct = gadd(ct, gmul(gen_2, gel(pk, i + 1)));
        }
        //ct = Fp_red(ct, gel(pk, 1));
    }
    //cout << GENtostr(ct) << endl;
    ct = Fp_red(ct, gel(pk, 1));
    return ct;
}

GEN encrypt(GEN sk, GEN pk, GEN m){
    //printf("Checkpoint: Encrypt");
    //cout<<GENtostr(m)<<endl;
    m = decimal_to_binary_rev(m);
    long n = lg(m) - 1;
    GEN ct = cgetg(n + 1, t_VEC);
    for(int i = 0; i < n; i++)
        gel(ct, i + 1) = encrypt_bit(sk, pk, gel(m, i + 1));
    return ct;
}

GEN decrypt_bit(GEN sk, GEN ct){
    GEN pt = Fp_red(ct, sk);
    pt = Fp_red(pt, gen_2);
    return pt;
}

GEN decrypt(GEN sk, GEN ct){
    int n = lg(ct) - 1;
    GEN pt = cgetg(n + 1, t_VEC);
    for(int i = 0; i < n; i++)
        gel(pt, i + 1) = decrypt_bit(sk, gel(ct, i + 1));
    return pt;
}

GEN multiply_bit(GEN ct_1, GEN ct_2, GEN pk){
    GEN result = gmul(ct_1, ct_2);
    for(int i = 0; i < gamma; i++){
        result = Fp_red(result, gel(pk, gamma + tau - i));
    }
    return result;
}

GEN add_bit(GEN ct_1, GEN ct_2, GEN x_0){
    return gmod(gadd(ct_2, ct_1), x_0);
}

GEN addition_gate(GEN ct_1, GEN ct_2, GEN pk, int n){
    if (n == -1)
    	int n = lg(ct_1) - 1;
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

GEN multiplication_utility(GEN result, GEN ct, GEN pk, int i){
	int n = lg(ct) - 1;
	GEN temp = cgetg(n + i + 1, t_VEC);
	for(int j = 0; j < i; j++)
		gel(temp, j + 1) = gen_0;
	for(int j = 0; j < n; j++)
		gel(temp, j + 1 + i) = gel(ct, i);
	result = addition_gate(result, temp, pk, n + i);
	return result;
}

GEN multiplication_gate(GEN ct_1, GEN ct_2, GEN pk){
	int n = lg(ct_1) - 1, m = lg(ct_2) - 1;
	//cout << GENtostr(ct_1) << endl;
	//cout << GENtostr(ct_2) << endl;
	GEN result = cgetg(m + n + 1, t_VEC);
	GEN temp = cgetg(n + 1, t_VEC);
	for(int i = 0; i < m + n; i++)
		gel(result, i + 1) = gen_0;
	for(int j = 0; j < n; j++)
		gel(result, j + 1) = multiply_bit(gel(ct_1, j + 1), gel(ct_2, 1), pk);
	//cout << GENtostr(result) << endl;
	for(int i = 1; i < m; i++){
		for(int j = 0; j < n; j++)
			gel(temp, j + 1) = multiply_bit(gel(ct_1, j + 1), gel(ct_2, i + 1), pk);
		//cout << GENtostr(temp) << " " << i << endl;
		result = multiplication_utility(result, temp, pk, i);
	}
	return result;
}

int main() {
    pari_init(600000000, 2);
    pari_sp ltop_super;
    GEN sk, pk;
    ifstream sk_file("secret_key.txt");
    if(!sk_file) {
        //printf("Error: Secret Key!");
        ltop_super = avma;
        sk = generate_secret_key();
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
        GEN pk = key_gen(sk);
        ofstream pk_file("public_key.txt");
        for(int i = 0; i < tau + gamma; i++)
            pk_file << GENtostr(gel(pk, i + 1)) << endl;
        pk_file.close();
    }
    else{
        pk = cgetg(tau + 1 + gamma, t_VEC);
        for(int i = 0; i < tau + gamma; i++){
            gel(pk, i + 1) = gp_read_stream(pk_file);
            //cout << GENtostr(gel(pk, i + 1)) << endl;
        }
        fclose(pk_file);
    }
    //printf("Hello!");
    //GEN ct_1 = encrypt_bit(sk, pk, stoi(0));
    //GEN ct_2 = encrypt_bit(sk, pk, stoi(1));

    GEN ct_1 = encrypt(sk, pk, stoi(2));
    GEN ct_2 = encrypt(sk, pk, stoi(3));


    GEN pt = decrypt(sk, multiplication_gate(ct_1, ct_2, pk));
    //GEN pt = decrypt(sk, multiplication_gate(ct_1, ct_2, pk));

    //GEN pt = decrypt_bit(sk, multiply_bit(ct_1, ct_2, pk));
    //GEN pt = decrypt_bit(sk, add_bit(ct_2, multiply_bit(ct_2, ct_2, pk), pk));

    //GEN pt = multiplication_gate(decimal_to_binary_rev(stoi(6)), decimal_to_binary_rev(stoi(4)), pk);
    
    cout << GENtostr(pt) << endl;
    /*for(int i = 0; i < lg(pt) - 1; i++)
        cout << GENtostr(gel(pt, i + 1)) << endl;
    */
    pari_close();

    return 0;
}