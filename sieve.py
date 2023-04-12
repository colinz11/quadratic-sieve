from math import pow, sqrt, exp, log, ceil


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


#returns a^n mod p FINISHED
def power(a, n, p):
    ret = 1
    while n > 0:
        if n % 2 == 1:
            ret = (ret * a) % p
        a = (a * a) % p
        n //= 2
    return ret


#get legendre number, we want when result is 1 FINISHED
def legendre(a, p):
    return power(a, (p - 1) // 2, p)

#generate primes up until B and use euler’s criterion to determine whether N is a quadratic residue mod p FINISHED
def find_factor_base(N, B):
    if B < 2:
        return []
    
    is_prime = [True] * (B+1)
    p = 2

    factor_base = []
    while p < B + 1: 
        if is_prime[p]:
            if legendre(N, p) == 1:
                factor_base.append(p)
            for j in range(2 * p, B + 1, p):
                is_prime[j] = False
        p += 1
    return factor_base
   
#picks optimal bound B FINISHED
def get_smoothness_bound(N):
    epsilon = 0.23
    B = exp((0.5 + epsilon) * sqrt(log(N)*log(log(N))))
    return int(B)

#does (x + n)^2 - N where x = sqrt(N) to find B-smooth numbers. Then use Tonelli and Shanks to compute resiudes for each prime in factor base
def find_smooth(factor_base, N, interval, tolerance=1):
    root = ceil(sqrt(N))
    sieve_static = [int(pow(root + i, 2) - N) for i in range(interval)]
 
    sieve_list = sieve_static.copy()

    prime_factors = {x : [] for x in sieve_static}
    
    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i, len(sieve_static), 2): # found the 1st even term, now every other term will also be even
            while sieve_list[j] % 2 == 0: #account for powers of 2
                sieve_list[j] //= 2
                prime_factors[sieve_static[j]].append(2)

    for p in factor_base[1:]: #skip 2
        residues = tonelli_shanks(N, p) 
        for r in residues: # two solutions r and p-r
            for i in range((r-root) % p, len(sieve_static), p): # Now every pth term will also be divisible
                while sieve_list[i] % p == 0:
                    sieve_list[i] //= p #divide by prime in factor base
                    prime_factors[sieve_static[i]].append(p)

    x = [] 
    smooth_nums = []
    factors = {}
    indices = [] 

    

    for i in range(len(sieve_static)):
        if len(smooth_nums) >= len(factor_base) + tolerance: 
            break
        if sieve_list[i] == 1: #found B smooth
            smooth_nums.append(sieve_static[i])
            factors[sieve_static[i]] = prime_factors[sieve_static[i]]
            x.append(i+root) #x+n
            indices.append(i) #n
    return smooth_nums, factors, x, indices

#builds matrix of exponents of prime factors of smooth numbers mod 2
def build_matrix(factor_base, smooth_nums, factors):
    matrix = []

    for n in smooth_nums:
        exp_vector = [0]*(len(factor_base))
        prime_factors = factors[n]

        for i in range(len(factor_base)):
            if factor_base[i] in prime_factors:
                exp_vector[i] = (exp_vector[i] + prime_factors.count(factor_base[i])) % 2 #get exponenet vector of prime factors mod 2
        matrix.append(exp_vector)
    print(matrix)
    return matrix

def transpose(matrix):
    return [[matrix[i][j] for i in range(len(matrix))] for j in range(len(matrix[0]))]

def gauss_elimination(matrix):
    pass


#newtons method for sqrt
def newton_sqrt(n): 
    return sqrt(n)



 

#is n prime with some probability 
def miller_rabin(n): 
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    pass

#tonelli algorithim copied from internet, two solutionrs r and p-r
def tonelli_shanks(n, p): 
    q = p - 1
    s = 0

    while q % 2 == 0:
        q //= 2
        s += 1

    if s == 1:
        r = power(n, (p + 1) // 4, p)
        return r, p-r
    
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break

    c = power(z, q, p)
    r = power(n, (q + 1) // 2, p)
    t = power(n, q, p)
    m = s
    t2 = 0

    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = power(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r, p-r

def main():
    N = 227179
    B = get_smoothness_bound(N)
    
    print(B)
    factor_base = find_factor_base(N, B)
    print(factor_base)
    smooth_nums, factors, x, indices = find_smooth(factor_base, N, B**3)
    if len(smooth_nums) < len(factor_base):
        print('No smooth numbers')
        return
    build_matrix(factor_base, smooth_nums, factors)


if __name__ == "__main__":
    main()
