from math import sqrt, exp, log, ceil
from itertools import chain
import time

#1. generate B-smooth factor base
#2. Find B-smooth relations using sieving and Tonelli-Shanks
#3. Build exponent matrix mod 2 from relations
#4. Solve matri for null space finding perfect squares
#5. Solve the congruence of squares to obtain factors 


#euclidian gcd algorithim
def gcd(a,b):
    if b > a:
        a, b = b, a

    while b != 0:
        a %= b
        if a == 0:
            return b
        b %= a
    return 1

#compute a^n
def int_pow(a, n):
    result = 1
    for i in range(n):
        result *= a
    return result



#returns a^n mod p
def power(a, n, p):
    ret = 1
    while n > 0:
        if n % 2 == 1:
            ret = (ret * a) % p
        a = (a * a) % p
        n //= 2
    return ret


#get legendre number, we want when result is 1
def legendre(a, p):
    return power(a, (p - 1) // 2, p)


#generate primes up until B and use euler’s criterion to determine whether N is a quadratic residue mod p
def find_factor_base(N, B):
    if B < 2:
        return []
    
    is_prime = [True] * (B+1)
    p = 2

    factor_base = []
    while p < B + 1: 
        if is_prime[p]:
            if legendre(N,p) == 1:
                    factor_base.append(p)
            for j in range(2 * p, B + 1, p):
                is_prime[j] = False
        p += 1
    return factor_base
   
#picks optimal bound B
def get_smoothness_bound(N):
    epsilon = 0.1
    e = 2.71
    B = e ** ((0.5 + epsilon) * sqrt(ceil(log(N)*log(log(N)))))
    return ceil(B)

#does (x + n)^2 - N where x = sqrt(N) to find B-smooth numbers.
def find_smooth(factor_base, N, interval, tolerance=3):
    root = ceil(sqrt(N))
    x = [] 
    smooth_nums = []
    factors = {}
    for i in range(interval):
        sieve_num = int(int_pow(root + i, 2) - N)
        sieve_static = [sieve_num][0]
        prime_factors = []
        for p in factor_base: 
            while sieve_num % p == 0:#p is a prime factor
                sieve_num //= p #divide by prime in factor base
                prime_factors.append(p)

        
        if sieve_num == 1:
            smooth_nums.append(sieve_static)
            factors[sieve_static] = prime_factors
            x.append(i+root) #x+n
            #if (len(smooth_nums) % 10 == 0):
            #   print(str(time.ctime()) + 'Smooth Number ' + str(len(smooth_nums))  + ':'+ str(sieve_static))
        if len(smooth_nums) >= (len(factor_base) * tolerance) + 1: 
              break
    return smooth_nums, x, factors

#builds matrix of exponents of prime factors of smooth numbers mod 2
def build_matrix(factor_base, smooth_nums, factors):
    matrix = [] #binary version
    matrix2 = [] #non-binary version

    for n in smooth_nums:
        exp_vector = [0]*(len(factor_base))
        exp_vector2 = [0]*(len(factor_base))
        prime_factors = factors[n]

        for i in range(len(factor_base)):
            if factor_base[i] in prime_factors:
                exp_vector[i] = (exp_vector[i] + prime_factors.count(factor_base[i])) % 2 #get exponenet vector of prime factors mod 2
                exp_vector2[i] = prime_factors.count(factor_base[i])
        matrix.append(exp_vector)
        matrix2.append(exp_vector2)
    return matrix, matrix2

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
def solve(N, solution_vector, smooth_nums, x, factor_base, matrix2):
    nums = [smooth_nums[i] for i in solution_vector]
    x_nums = [x[i] for i in solution_vector]

    a = 1
    b = 1
    smooth_vector = [0 for i in factor_base] #this is an exponent vector


    for i in solution_vector:
        for j in range(len(factor_base)):
            smooth_vector[j] += matrix2[i][j]

    for i in range(len(smooth_vector)):
        f = 1
        for j in range(smooth_vector[i] >> 1):
            f = f * (factor_base[i])

        a *= f

    for n in x_nums:
        b *= n


    
    return gcd(a+b, N)    

def main():
    #N = 1009 * 191161
    #N = 4201 * 8011
    #N = 484459 * 191161
    #N = 4201 * 484459
    #N = 1000033 * 1000003
    #N = 10000019 * 10000169
    #N = 100000007 * 10000169
    #N = 100000007 * 100000049
    #N = 16921456439215439701
    N = 46839566299936919234246726809
    #N = 6172835808641975203638304919691358469663

    B = get_smoothness_bound(N)

    print('Bound: ' + str(B))
    
    factor_base = find_factor_base(N, B)
    print("num primes in factor base is " + str(len(factor_base)))

    
    print('Finding smooth numbers ...')
    smooth_nums, x, factors = find_smooth(factor_base, N, B**3, tolerance=1.01)

    if len(smooth_nums) < len(factor_base):
        print('No smooth numbers')
        return
    
    print('Found smooth numbers')
    #below, matrix is all mod 2 entries, while matrix2 has any positive integer for each entry
    matrix, matrix2 = build_matrix(factor_base, smooth_nums, factors) 

    print('Running guass elimination ... ')
    solution_rows, marks, matrix = gauss_elimination(matrix)



    for K in range(len(solution_rows)):
        print('Finding factors ...')
        solution_vector = solve_row(matrix, marks, solution_rows, K)
        factor = solve(N, solution_vector, smooth_nums, x, factor_base, matrix2)
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
