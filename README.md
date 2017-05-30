(Leveled) Homomorphic C++ implementation of DGHV scheme using GNU MP library and Pari C library.
The implementation supports addition and multiplication gates on integers.

For the parameters taken in the implementation, the SwHE scheme supports approximately 350 multiplications.

I've used rand(), seeded with system time to generate random numbers, which is cryptographically insecure.
But, this is not an industrial grade implementation.

For the current parameters, here's an example of recrypt (ignore the degree after recrypt, need to fix that):
![](./Recrypt_Example.png)
