from math import sqrt, log, ceil
from itertools import chain
import time
import random

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

#compute a^n
def int_pow(a, n):
    result = 1
    for i in range(n):
        result *= a
    return result



#returns a^n mod p FINISHED
def exp_mod(a, n, p):
    ret = 1
    while n > 0:
        if n % 2 == 1:
            ret = (ret * a) % p
        a = (a * a) % p
        n //= 2
    return ret


#get legendre number, we want when result is 1 FINISHED
def legendre_num(a, p):
    return exp_mod(a, (p - 1) // 2, p)


#generate primes up until B and use euler’s criterion to determine whether N is a quadratic residue mod p FINISHED
def find_factor_base(N, B):
    if B < 2:
        return []
    
    prime_bool = [True for i in range(B + 1)]
    p = 2

    factor_base = []
    while p < (B + 1): 
        if prime_bool[p]:
            if legendre_num(N,p) == 1:
                    factor_base.append(p)

            #mark off all multiples of p as being composite
            for j in range(2 * p, B + 1, p):
                prime_bool[j] = False
        p += 1
    return factor_base
   
#picks optimal bound B FINISHED
def set_B(N, epsilon=0.1):
    e = 2.71828183
    B = e ** ((0.5 + epsilon) * sqrt(ceil(log(N)*log(log(N)))))
    return ceil(B)

#does (x + n)^2 - N where x = sqrt(N) to find B-smooth numbers. Then use Tonelli and Shanks to compute resiudes for each prime in factor base
def find_smooth_numbers(factor_base, N, interval, extra_rows=1):
    root = ceil(sqrt(N))
    x = [] 
    smooth_nums = []
    factors = {}
    for i in range(interval):
        
        #sieve_num = int(int_pow(root + i, 2) - N)
        sieve_num = int((root + i) * (root + i) - N)
        if sieve_num > N:
            sieve_num = sieve_num % N
        if sieve_num == 0:
            continue
    
  
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
            if (len(smooth_nums) % 10 == 0):
                print('Smooth Number ' + str(len(smooth_nums))  + ':'+ str(sieve_static))
        if len(smooth_nums) >= (len(factor_base) + extra_rows): 
              break
    return smooth_nums, x, factors



#builds matrix of exponents of prime factors of smooth numbers mod 2
def build_matrix(factor_base, smooth_nums, factors):
    matrix = [] #binary version
    matrix_nb = [] #non-binary version

    for n in smooth_nums:
        exp_vector = [0]*(len(factor_base))
        exp_vector2 = [0]*(len(factor_base))
        prime_factors = factors[n]

        for i in range(len(factor_base)):
            if factor_base[i] in prime_factors:
                exp_vector[i] = (exp_vector[i] + prime_factors.count(factor_base[i])) % 2 #get exponenet vector of prime factors mod 2
                exp_vector2[i] = prime_factors.count(factor_base[i])
        matrix.append(exp_vector)
        matrix_nb.append(exp_vector2)
    return matrix, matrix_nb


def matrix_transpose(matrix):
    original_num_rows = len(matrix)
    original_num_cols = len(matrix[0])
    new_matrix = [[0 for i in range(original_num_rows)] for j in range(original_num_cols)]
    for i in range(original_num_cols):
        for j in range(original_num_rows):
            new_matrix[i][j] = matrix[j][i]
    return new_matrix

def find_left_null_space(matrix):
    num_rows = len(matrix)
    num_cols = len(matrix[0])
    solution_combos = [[0 for i in range(num_rows)] for j in range(num_rows)]
    rows = [i for i in range(num_rows)]
    for i in range(num_rows):
        solution_combos[i][i] = 1
    for piv_row_index in range(num_cols):
        if (piv_row_index % 100 == 0):
            print("piv row index " + str(piv_row_index))
        piv_row = matrix[piv_row_index]
        #need to ensure that piv_row[piv_row_index] == 1
        if (piv_row[piv_row_index] != 1):
            new_row_index = piv_row_index + 1
            while (new_row_index < num_rows):
                if matrix[new_row_index][piv_row_index] == 1:
                    #we can swap this with the pivot row
                    original_copy = [k for k in piv_row]
                    new_copy = [k for k in matrix[new_row_index]]
                    matrix[piv_row_index] = new_copy
                    matrix[new_row_index] = original_copy
                    piv_row = new_copy

                    row_copy = rows[piv_row_index]
                    rows[piv_row_index] = rows[new_row_index]
                    rows[new_row_index] = row_copy
                    break
                new_row_index += 1
        if piv_row[piv_row_index] != 1:
            continue
        #other_row_indices = [i for i in range(piv_row_index)] + [i for i in range(piv_row_index +1, num_rows)]     
        other_row_indices =  [i for i in range(piv_row_index +1, num_rows)]   
        #print(other_row_indices)       
        for other_row_index in other_row_indices:
   
            other_row = matrix[other_row_index]
            if (other_row[piv_row_index] == 1):
                #add piv_row (mod 2) to get a 0 in the pivot column
                matrix[other_row_index] = [(other_row[i] + piv_row[i]) % 2 for i in range(num_cols)]

                row1 = rows[piv_row_index]
                row2 = rows[other_row_index]
                solution_combos[row2] = [(solution_combos[row1][i] + solution_combos[row2][i]) % 2 for i in range(num_rows)]
    
    
    return [solution_combos[rows[i]] for i in range(num_cols, num_rows)]      




#gives list of null spaces
def gauss_elimination(matrix):
    matrix = matrix_transpose(matrix)
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
            
  
    matrix = matrix_transpose(matrix)
    
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
def find_squared_congruence(N, solution_vector, smooth_nums, x, factor_base, matrix_nb):
    solution_rows = []
    for i in range(len(solution_vector)):
        if solution_vector[i] == 1:
            solution_rows.append(i)
    nums = [smooth_nums[i] for i in solution_rows]
    x_nums = [x[i] for i in solution_rows]
    

    #print(nums)
    #print(x_nums)
    a = 1
    b = 1
    smooth_vector = [0 for i in factor_base] #this is an exponent vector


    for i in solution_rows:
        for j in range(len(factor_base)):
            smooth_vector[j] += matrix_nb[i][j]
  
    for i in range(len(smooth_vector)):
        f = 1
        for j in range(smooth_vector[i] >> 1):
            f = f * (factor_base[i])
        a *= f

    for n in x_nums:
        b *= n
    return (a, b)   



# run Miller rabin several times with different values of a
def check_prime(n, iter=10):
    if (pow(2, n-1, n) != 1):
        return False
    m = 1
    beta = n -1
    k = 0
    while (beta % 2 == 0):
        beta = beta >> 1
        k += 1
    m = (n-1) >> k 
    for i in range(iter):
        a = random.randint(2, n-2)
        b_0 = exp_mod(a, m, n)
        if (b_0 == 1) or (b_0 == -1) or (b_0 == n-1):
            #n is probably prime
            continue
        b_j = b_0
        j = 1
        while (j) < k:
            b_j = (b_j * b_j) % n
            if (b_j == 1):
                return False
            elif (b_j == -1) or (b_j == n-1):
                #n is probably prime
                break
            j += 1
    return True

#tonelli algorithim copied from internet, two solutionrs r and p-r
def tonelli_shanks(n, p): 
    q = p - 1
    s = 0

    while q % 2 == 0:
        q //= 2
        s += 1

    if s == 1:
        r = exp_mod(n, (p + 1) // 4, p)
        return r, p-r
    
    for z in range(2, p):
        if p - 1 == legendre_num(z, p):
            break

    c = exp_mod(z, q, p)
    r = exp_mod(n, (q + 1) // 2, p)
    t = exp_mod(n, q, p)
    m = s
    t2 = 0

    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = exp_mod(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r, p-r

def main(N=21, epsilon=0.1):
    start_time_ms = int(time.time() * 1000)
    #N = 55587
    #N = 1009 * 191161
    #N = 37 * 1009
    #N = 484459 * 191161
    #N = 4201 * 484459
    #N = 1000033 * 1000003
    #N = 10000019 * 10000169
    #N = 100000007 * 10000169
    #N = 100000007 * 100000049
    #N = 16921456439215439701
    #N = 100000007 * 10000169312

    N = 46839566299936919234246726809
    #N = 6172835808641975203638304919691358469663
    if check_prime(N):
        print("N is very likely prime")
        return

    B = set_B(N, epsilon)

    print('Bound: ' + str(B))
    
    factor_base = find_factor_base(N, B)
    #print("num primes in factor base is " + str(len(factor_base)))
    #print("factor base is " + str(factor_base))

    time_0 = int(time.time() * 1000)
    print('Finding smooth numbers ...')
    smooth_nums, x, factors = find_smooth_numbers(factor_base, N, B**3, extra_rows=20)
    #print("smooth_nums is " + str(smooth_nums))
    #print("x is " + str(x))
    found_smooth_time_ms = int(time.time() * 1000) - time_0

    if len(smooth_nums) < len(factor_base):
        print('Not enough B-smooth numbers were found')
        return
    
    print('Found smooth numbers')
    #below, matrix is all mod 2 entries, while matrix_nb has any positive integer for each entry
    matrix, matrix_nb = build_matrix(factor_base, smooth_nums, factors) 


    time_1 = int(time.time() * 1000)
    print('Running Gaussian elimination ... ')
    #solution_rows, marks, matrix = gauss_elimination(matrix)
    

    solution_combos = find_left_null_space(matrix)
    elim_time = int(int(time.time() * 1000) - time_1)
    num_solutions = len(solution_combos)



    #for K in range(len(solution_rows)):
    for K in range(num_solutions):
        print('Finding factors ...')
        #solution_vector = solve_row(matrix, marks, solution_rows, K)
        #[a, b] = find_squared_congruence(N, solution_vector, smooth_nums, x, factor_base, matrix_nb)
        [a, b] = find_squared_congruence(N, solution_combos[K], smooth_nums, x, factor_base, matrix_nb)
        factor = abs(gcd(a-b, N))

        if factor == 1 or factor == N:
            print("Trivial Factor")
            continue
        else:
            print("N = " + str(factor) + " * " + str(int(N/factor)))
            total_time = int(time.time() * 1000) - start_time_ms
            return [found_smooth_time_ms, elim_time, total_time]
        
    #if the code reaches this point, no nontrivial factors have been found
    print("Could not find a nontrivial factor")

if __name__ == "__main__":
    main()
