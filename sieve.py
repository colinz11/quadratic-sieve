from math import sqrt, log, ceil
from itertools import chain
import time
import random

#1. generate B-smooth factor base
#2. Find B-smooth relations using sieving
#3. Build exponent matrix mod 2 from relations
#4. Solve matri for null space finding perfect squares
#5. Solve the congruence of squares to obtain factors 


#euclidian gcd algorithim
def gcd(a,b):
    if b > a: #make it so a > b
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
def exp_mod(a, n, p):
    ret = 1
    while n > 0:
        if n % 2 == 1:
            ret = (ret * a) % p
        a = (a * a) % p
        n //= 2
    return ret

#get legendre number, we want when result is 1
def legendre_num(a, p):
    return exp_mod(a, (p - 1) // 2, p)

#generate primes up until B and use euler’s criterion to determine whether N is a quadratic residue mod p
def find_factor_base(N, B):
    if B < 2:
        return []
    
    prime_bool = [True for i in range(B + 1)]
    p = 2 # start a first prime

    factor_base = []
    while p < (B + 1): 
        if prime_bool[p]:
            if legendre_num(N,p) == 1: #factor base only includes primes with legrendre num equal to 1
                    factor_base.append(p) 

            #mark off all multiples of p as being composite
            for j in range(2 * p, B + 1, p):
                prime_bool[j] = False
        p += 1
    return factor_base
   
#picks optimal bound B
def set_B(N, epsilon=0.1): 
    e = 2.71828183
    B = e ** ((0.5 + epsilon) * sqrt(ceil(log(N)*log(log(N)))))
    return ceil(B)

#does (x + n)^2 - N where x = sqrt(N) to find B-smooth numbers.
def find_smooth_numbers(factor_base, N, interval, extra_rows=1):
    root = ceil(sqrt(N))
    x = [] 
    smooth_nums = []
    factors = {}
    for i in range(interval):
        
        #f_i = int(int_pow(root + i, 2) - N)
        f_i = int((root + i) * (root + i) - N)
        if f_i > N:
            f_i = f_i % N

        if f_i == 0:
            #this would create an infinite loop a few lines later
            continue
    
        f_i_original = f_i
        prime_factors = []
        for p in factor_base: 
            while f_i % p == 0:#p is a prime factor
                f_i //= p #divide by prime in factor base
                prime_factors.append(p)

        
        if f_i == 1:
            smooth_nums.append(f_i_original)
            factors[f_i_original] = prime_factors
            x.append(i+root) #x+n
            if (len(smooth_nums) % 10 == 0):
                print('Smooth Number ' + str(len(smooth_nums))  + ':'+ str(f_i_original))
        if len(smooth_nums) >= (len(factor_base) + extra_rows): 
              break
    return smooth_nums, x, factors

def lower_register(orig, p):
    while (orig % p) == 0:
        orig //= p
    return orig

def lower_register_log(orig, log_p):
    return orig - log_p

def check_smooth(f, x, N, factor_base, smooth_nums, x_list, factors):
    if f < 1:
        return
    if f == 1:
        #f(x) is smooth, where f(x) = (x)^2 - N
        x_list.append(x)
        f_orig = x * x - N
        smooth_nums.append(f_orig)
        if (len(smooth_nums) % 20 == 0):
            print("smooth nums: " + str(len(smooth_nums)))
        #now find factors of f_orig
        prime_factors = []
        for p in factor_base:
            while f_orig % p == 0:
                f_orig //= p
                prime_factors.append(p)
        factors[x * x - N] = prime_factors

def check_smooth_log(f, x, N, factor_base, smooth_nums, x_list, factors, threshold):
    if f <= 0:
        return
    elif f < threshold:

        #f(x) is smooth, where f(x) = (x)^2 - N
        '''
        x_list.append(x)
        f_orig = x * x - N
        smooth_nums.append(f_orig)
        if (len(smooth_nums) % 20 == 0):
            print("smooth nums: " + str(len(smooth_nums)))
        #now find factors of f_orig
        prime_factors = []
        for p in factor_base:
            while f_orig % p == 0:
                f_orig //= p
                prime_factors.append(p)
        factors[x * x - N] = prime_factors
        '''
        prime_factors = []
        f_orig = x * x - N
        for p in factor_base:
            while f_orig % p == 0:
                f_orig //= p
                prime_factors.append(p)
        if (f_orig > 1):
            return
        else:
            if (len(smooth_nums) % 20 == 0):
                print(str(time.ctime()) + ": " + str(len(smooth_nums)) + " smooths have been found")
            x_list.append(x)
            smooth_nums.append(x * x - N)
            factors[x * x - N] = prime_factors
        






def find_smooth_numbers_TS(factor_base, N, interval, extra_rows=1):
    root = ceil(sqrt(N))
    x_list = [] 
    smooth_nums = []
    factors = {}
    num_needed = len(factor_base) + extra_rows
    print("number of smooth numbers needed is " + str(num_needed))

    solutions_mod_p = [TS(i, N) for i in factor_base]

    x_upper_lim = int((sqrt(2) - 1)*sqrt(N) - 1)
   

    length = min(max(factor_base) * 10, x_upper_lim)
    log_p = [log(p) for p in factor_base]
    threshold = log(max(factor_base) + 1)
    
    end=-1
    while len(smooth_nums) < num_needed:
        start=end+1
        end = start+length

        #f_x = [(root + start + i) * (root + start + i) - N for i in range(length)]
        f_x = [log((root + start + i) * (root + start + i) - N) for i in range(length)]
        for i in range(len(factor_base)):
            p_i = factor_base[i]
            (x1, x2) = solutions_mod_p[i]
            #iterate through the f_x array dividing out p_i 
            #first_val_1 = (ceil((root + start) / p_i)  * p_i) - start - root + x1 #first val s.t. f_x[first_val_1] divisible by p_i
            first_val_1 = (0 - ((root + start) // (0 - p_i))  * p_i) - start - root + x1
            v1 = root + start
            v2 = 0 -(v1 // (0 - p_i))
            v3 = ceil(v2)
            v4 = v3 * p_i
            v5 = v4 - start - root + x1
            #first_val_2 = (ceil((root + start) / p_i)  * p_i) - start - root + x2
            first_val_2 = (0 - ((root + start) // (0 - p_i))  * p_i) - start - root + x2
            if first_val_1 >= p_i:
                first_val_1 -= p_i
            if first_val_2 >= p_i:
                first_val_2 -= p_i


            
            for x_val in range(first_val_1, length, p_i):
                #divide out p
                #f_x[x_val] = lower_register(f_x[x_val], p_i)
                f_x[x_val] = lower_register_log(f_x[x_val], log_p[i])
                #check that x_val+start is correct argument
                #check_smooth(f_x[x_val], x_val+start+root, N, factor_base, smooth_nums, x_list, factors)
                check_smooth_log(f_x[x_val], x_val+start+root, N, factor_base, smooth_nums, x_list, factors, threshold)

            for x_val in range(first_val_2, length, p_i):
                #divide out p
                #f_x[x_val] = lower_register(f_x[x_val], p_i)
                f_x[x_val] = lower_register_log(f_x[x_val], log_p[i])
                #check_smooth(f_x[x_val], x_val+start+root, N, factor_base, smooth_nums, x_list, factors)
                check_smooth_log(f_x[x_val], x_val+start+root, N, factor_base, smooth_nums, x_list, factors, threshold)


    return smooth_nums, x_list, factors




#builds matrix of exponents of prime factors of smooth numbers mod 2
def make_matrix(factor_base, smooth_nums, factors):
    matrix = [] #binary version
    matrix_nb = [] #non-binary version

    for n in smooth_nums:
        exp_vector = [0 for i in range(len(factor_base))]
        exp_vector2 = [0 for i in range(len(factor_base))]
        prime_factors = factors[n]

        for i in range(len(factor_base)):
            if factor_base[i] in prime_factors:
                exp_vector[i] = (exp_vector[i] + prime_factors.count(factor_base[i])) % 2 #get exponenet vector of prime factors mod 2
                exp_vector2[i] = prime_factors.count(factor_base[i])
        matrix.append(exp_vector)
        matrix_nb.append(exp_vector2)
    return matrix, matrix_nb

#tranposes a matrix
def matrix_transpose(matrix):
    original_num_rows = len(matrix)
    original_num_cols = len(matrix[0])
    new_matrix = [[0 for i in range(original_num_rows)] for j in range(original_num_cols)]
    for i in range(original_num_cols):
        for j in range(original_num_rows):
            new_matrix[i][j] = matrix[j][i]
    return new_matrix

#doess row reductions to find the left null space
def find_left_null_vectors(matrix):
    num_rows = len(matrix)
    num_cols = len(matrix[0])
    combos = [[0 for i in range(num_rows)] for j in range(num_rows)]
    rows = [i for i in range(num_rows)]
    for i in range(num_rows):
        combos[i][i] = 1
    for piv_row_index in range(num_cols):
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
                combos[row2] = [(combos[row1][i] + combos[row2][i]) % 2 for i in range(num_rows)]
    
    return [combos[rows[i]] for i in range(num_cols, num_rows)]      

#caulates a^2 = b^2 mod n and does gcd(a+b, n)
def find_squared_congruence(null_vector, x, factor_base, matrix_nb):
    combo_rows = []
    for i in range(len(null_vector)):
        if null_vector[i] == 1:
            combo_rows.append(i)
    x_nums = [x[i] for i in combo_rows]
    
    a = 1
    b = 1
    smooth_vector = [0 for i in factor_base] #this is an exponent vector

    for i in combo_rows:
        for j in range(len(factor_base)):
            smooth_vector[j] += matrix_nb[i][j] #gets the exponent power for each prime
  
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


#copied from internet
def TS(p, n):

    if legendre_num(n, p) != 1:
        print("not a square (mod p)")
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        r = pow(n, (p + 1) // 4, p)
        r2 = (p - r) % p
        return (r, r2)
    z = 2
    for z_val in range(2, p):
        z = z_val
        if p - 1 == legendre_num(z_val, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return (r, (p-r)%p)

def factor(N=21, epsilon=0.001):
    start_time_ms = int(time.time() * 1000)

    if check_prime(N):
        print("N is very likely prime")
        return

    B = set_B(N, epsilon)

    print('Bound: ' + str(B))
    
    factor_base = find_factor_base(N, B)

    time_0 = int(time.time() * 1000)
    print('Finding smooth numbers ...')
    extra_rows = int(20 + log(N))
    if (N < ( 10 ** 15 )):
        smooth_nums, x, factors = find_smooth_numbers(factor_base, N, B**3, extra_rows)
    else:
        smooth_nums, x, factors = find_smooth_numbers_TS(factor_base, N, B**3, extra_rows)
        if smooth_nums == None:
            smooth_nums, x, factors = find_smooth_numbers(factor_base, N, B**3, extra_rows)
    found_smooth_time_ms = (time.time() * 1000) - time_0

    if len(smooth_nums) < len(factor_base):
        print('Not enough B-smooth numbers were found')
        return
    else:
        print('Found smooth numbers in ' + str(found_smooth_time_ms / 1000) + "s")

    #below, matrix is all mod 2 entries, while matrix_nb has any positive integer for each entry
    matrix, matrix_nb = make_matrix(factor_base, smooth_nums, factors) 

    time_1 = int(time.time() * 1000)
    print('Running Gaussian elimination ... ')
    
    null_combos = find_left_null_vectors(matrix)
    elim_time = int(int(time.time() * 1000) - time_1)
    print("elim time: " + str(elim_time / 1000) + "s")
    num_solutions = len(null_combos)

    for s in range(num_solutions):
        [a, b] = find_squared_congruence(null_combos[s], x, factor_base, matrix_nb)
        possible_factor = abs(gcd(a-b, N))

        if possible_factor == 1 or possible_factor == N:
            print("Trivial Factor")
            continue
        else:
            print("N = " + str(possible_factor) + " * " + str(int(N/possible_factor)))
            total_time = int(time.time() * 1000) - start_time_ms
            print("total time: " + str(total_time) + "ms, finding smooths time: " + str(found_smooth_time_ms) +
                   "ms, Gaussian elim time: " + str(elim_time) + "ms")
            return [found_smooth_time_ms, elim_time, total_time]
        
    #if the code reaches this point, no nontrivial factors have been found
    #retry with higher epsilon
    print("Could not find a nontrivial factor, retrying to factor")
    return factor(N, epsilon * 10)

def main():
    #N= 21
    #N = 4201 * 8011
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
    N = 46839566299936919234246726809
    #N = 1000000000100011 * 1084051191974761
    #N = 6172835808641975203638304919691358469663
    #N = 3744843080529615909019181510330554205500926021947
    return factor(N, epsilon=0.005)
if __name__ == "__main__":
    print(main())

