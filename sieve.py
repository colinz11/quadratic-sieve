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


#returns a^n mod p 
def power(a, n, p):
    ret = 1
    while n > 0:
        if n % 2 == 1:
            ret = (ret * a) % p
        a = (a * a) % p
        n //= 2
    return ret

#generate primes up until B and use euler’s criterion to determine whether N is a quadratic residue mod p FINISHED
def find_factor_base(N, B):
    if B < 2:
        return []
    
    is_prime = [True] * (B+1)
    p = 2

    factor_base = []
    while p * p < B: 
        if is_prime[p]:
            if legendre(N, p) == 1:
                factor_base.append(p)
            for j in range(p * p, B + 1, p):
                is_prime[j] = False
        p += 1
    return factor_base
   
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
    return power(a, (p - 1) // 2, p)
 

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
    N = 124871214124141680743806720389
    print(find_factor_base(N, get_smoothness_bound(N)))

if __name__ == "__main__":
    main()
