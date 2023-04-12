from math import pow, sqrt, exp, log, ceil
from itertools import chain

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
            if legendre(N, p) == 1: #finds our factor base
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

   
    
    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i, len(sieve_static), 2): #every other term is now even
            while sieve_list[j] % 2 == 0: #2 is a prime factor
                sieve_list[j] //= 2
               

    for p in factor_base[1:]: #skip 2
        residues = tonelli_shanks(N, p) 
        for r in residues: # two solutions r and p-r
            for i in range((r-root) % p, len(sieve_static), p): #every pth term
                while sieve_list[i] % p == 0:#p is a prime factor
                    sieve_list[i] //= p #divide by prime in factor base


    x = [] 
    smooth_nums = []

    

    for i in range(len(sieve_static)):
        if len(smooth_nums) >= len(factor_base) + tolerance: 
            break
        if sieve_list[i] == 1: #found B smooth
            smooth_nums.append(sieve_static[i])

            x.append(i+root) #x+n
    return smooth_nums, x

#builds matrix of exponents of prime factors of smooth numbers mod 2
def build_matrix(factor_base, smooth_nums, factors):
    matrix = []

    for n in smooth_nums:
        exp_vector = [0]*(len(factor_base))
        prime_factors = factor(n, factor_base)

        for i in range(len(factor_base)):
            if factor_base[i] in prime_factors:
                exp_vector[i] = (exp_vector[i] + prime_factors.count(factor_base[i])) % 2 #get exponenet vector of prime factors mod 2
        matrix.append(exp_vector)
    return matrix

def factor(n, factor_base):#trial division from factor base
    factors = []

    for p in factor_base:
        while n % p == 0:
            factors.append(p)
            n //= p
    return factors


def transpose(matrix):
    return [[matrix[i][j] for i in range(len(matrix))] for j in range(len(matrix[0]))]

#gives list of null spaces
def gauss_elimination(matrix):
    matrix = transpose(matrix)
    marks = [False] * len(matrix[0])
    
    for i in range(len(matrix)): 
        row = matrix[i]
        for num in row: #search for pivot
            if num == 1:
                j = row.index(num) #column index
                marks[j] = True
                
                for k in chain(range(0, i), range(i + 1, len(matrix))): #search for other 1s
                    if matrix[k][j] == 1:
                        for i in range(len(matrix[k])):
                            matrix[k][i] = (matrix[k][i] + row[i]) % 2
                break
            
  
    matrix = transpose(matrix)
    
    solution_rows = []
    for i in range(len(marks)): #find free columns
        if not marks[i]:
            free_row = [matrix[i], i]
            solution_rows.append(free_row)
    return solution_rows, marks, matrix

#gets the correct indices from the null space
def solve_row(matrix, marks, solution_rows, K):
    
    row = solution_rows[K][0]
    indices = []

    solution_vector = []

    for i in range(len(row)):
        if row[i] == 1:
            indices.append(i)
    
    for r in range(len(matrix)):
        for c in indices:
            if matrix[r][c] == 1 and marks[r]:
                solution_vector.append(r)
    solution_vector.append(solution_rows[K][1])
    return solution_vector

#caulates a^2 = b^2 mod n and does gcd(a+b, n)
def solve(N, solution_vector, smooth_nums, x):
    nums = [smooth_nums[i] for i in solution_vector]
    x_nums = [x[i] for i in solution_vector]


    a = 1
    b = 1

    for n in nums:
        a *= n

    for n in x_nums:
        b *= n
    

    a = sqrt(a)


    return gcd(a+b, N)    


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
    N = 4423 * 6079
    B = get_smoothness_bound(N)

    print('Bound: ' + str(B))
    
    factor_base = find_factor_base(N, B)

    
    print('Finding smooth numbers ...')

    smooth_nums, x = find_smooth(factor_base, N, B**3)

    
    if len(smooth_nums) < len(factor_base):
        print('No smooth numbers')
        return
    
    print('Found smooth numbers')

    matrix = build_matrix(factor_base, smooth_nums)


    print('Running guass elimination ... ')
    solution_rows, marks, matrix = gauss_elimination(matrix)


    
    for K in range(len(solution_rows)):
        print('Finding factors ...')
        solution_vector = solve_row(matrix, marks, solution_rows, K)
        factor = solve(N, solution_vector, smooth_nums, x)
        if factor == 1 or factor == N:
            print(factor)
            print("Trivial factor")
            continue
        else:
            print(factor)
            print(N/factor)
            return

if __name__ == "__main__":
    main()
