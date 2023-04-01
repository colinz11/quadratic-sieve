from math import pow, sqrt


#1. generate B-smooth factor base
#2. Find B-smooth relations using sieving and Tonelli-Shanks
#3. Build exponent matrix mod 2 from relations
#4. Solve matri for null space finding perfect squares
#5. Solve the congruence of squares to obtain factors 

# https://en.wikipedia.org/wiki/Quadratic_sieve 

#euclidian gcd algorithim
def gcd(n): 
    pass

#newtons method for sqrt
def isqrt(n): 
    return sqrt(n)

def gauss_elimination(matrix):
    return matrix

def legendre(a, p):
    return 0 if a % p == 0 else pow(a, (p-1) // 2 , p) #what is pow

#is n prime with some probability 
def miller_rabin(n): 
    return False

#tonelli algorithim 
def tonelli(n, p): 
    pass

