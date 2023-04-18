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
    while b != 0:
        t = b
        b = a % b
        a = t

    return a


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
    epsilon = 0.5
    B =  exp((0.5 + epsilon) * sqrt(ceil(log(N)*log(log(N)))))
    return ceil(B)

#does (x + n)^2 - N where x = sqrt(N) to find B-smooth numbers. Then use Tonelli and Shanks to compute resiudes for each prime in factor base
def find_smooth(factor_base, N, interval, tolerance=1):
    root = ceil(sqrt(N))
    x = [] 
    smooth_nums = []
    factors = {}
    for i in range(interval):
        sieve_num = int(pow(root + i, 2) - N)
        sieve_static = [sieve_num][0]
        prime_factors = []
        for p in factor_base: 
            while sieve_num % p == 0: #p is a prime factor
                sieve_num //= p #divide by prime in factor base
                prime_factors.append(p)

        
        if sieve_num == 1:
            smooth_nums.append(sieve_static)
            factors[sieve_static] = prime_factors
            x.append(i+root) #x+n
            #print('Smooth Number: ' + str(sieve_static))
        
        if len(smooth_nums) >= len(factor_base) + tolerance: #stop when we find smooth numbers > # of primes
              break
    return smooth_nums, x, factors

#builds matrix of exponents of prime factors of smooth numbers mod 2
def build_matrix(factor_base, smooth_nums, factors):
    matrix = []

    for n in smooth_nums:
        exp_vector = [0]*(len(factor_base))
        prime_factors = factors[n] #factors for smooth number

        for i in range(len(factor_base)):
            if factor_base[i] in prime_factors:
                exp_vector[i] = (exp_vector[i] + prime_factors.count(factor_base[i])) % 2 #get exponenet vector of prime factors mod 2
        matrix.append(exp_vector)
    return matrix

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
        a *= n #calculating a^2
    for n in x_nums:
        b *= n

  
    a = newton_sqrt(a)
   
    return gcd(b+a, N)    

def congruent(a, b, n):
    return a % n == b % n

#newtons method for sqrt
def newton_sqrt(n): 
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

#is n prime with some probability 
def miller_rabin(n): 
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    pass

def main():
    N = 16921456439215439701
    
    B = get_smoothness_bound(N)

    print('Bound: ' + str(B))
    
    factor_base = find_factor_base(N, B)

    
    print('Finding smooth numbers ...')
    smooth_nums, x, factors = find_smooth(factor_base, N, B**3)

    if len(smooth_nums) < len(factor_base):
        print('Not enough smooth numbers')
        return
    
    print(f'Found {len(smooth_nums)} smooth numbers')

    print('Building matrix ...')
    matrix = build_matrix(factor_base, smooth_nums, factors)


    print('Running guass elimination ... ')
    solution_rows, marks, matrix = gauss_elimination(matrix)


    print('Finding factors ...')
    for K in range(len(solution_rows)):
        solution_vector = solve_row(matrix, marks, solution_rows, K)
        factor = solve(N, solution_vector, smooth_nums, x)
        if factor == 1 or factor == N:
            print("Trivial factor: " + str(factor))
        else:
            print(f"Found Factors: {factor}, {int(N/factor)}")
            return

if __name__ == "__main__":
    main()