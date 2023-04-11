from math import pow, sqrt, exp, log


#1. generate B-smooth factor base
#2. Find B-smooth relations using sieving and Tonelli-Shanks
#3. Build exponent matrix mod 2 from relations
#4. Solve matri for null space finding perfect squares
#5. Solve the congruence of squares to obtain factors 

# https://en.wikipedia.org/wiki/Quadratic_sieve 

#euclidian gcd algorithim FINISHED
def gcd(a,b):

    if b > a:
        a, b = b, a

    while b != 0:
        a %= b
        if a == 0:
            return b
        b %= a
    return 1

#generate all prime numbers up to n FINISHED
def generate_primes(n):
    if n < 2:
        return []
    
    is_prime = [True] * (n+1)
    i = 2

    while i * i < n: 
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = False
        i += 1

    return [p for p in range(2, n+1) if is_prime[p]]
 
#picks optimal bound B FINISHED
def get_smoothness_bound(N):
    epsilon = 0.1
    B = exp((0.5 + epsilon) * sqrt(log(N)*log(log(N))))
    return int(B)
    

#newtons method for sqrt
def isqrt(n): 
    return sqrt(n)

def gauss_elimination(matrix):
    return matrix

def legendre(a, p):
    return 0 if a % p == 0 else pow(a, (p-1) // 2 , p) #what is pow

#is n prime with some probability 
def miller_rabin(n): 
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    pass

#tonelli algorithim 
def tonelli(n, p): 
    pass

def main():
    print(get_smoothness_bound(1000000000000))

if __name__ == "__main__":
    main()
