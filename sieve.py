from math import sqrt, log, ceil
import time
import random

'''
Math 404 Cryptography project. 

A quadratic sieve algorithm to factor large numbers.

Authors: Keith Cressman, Aakash Kothapally, Colin Zhu
'''

def gcd(a,b):
    """Finds gcd of two numbers

    Args:
        a (int)
        b (int)

    Returns:
        int: Gcd of a and b
    """

    #make it so a > b
    if b > a: 
        a, b = b, a

    while b != 0:
        a %= b
        if a == 0:
            return b
        b %= a
    return 1

def exp_mod(a, n, p):
    """Calculates a^n mod p

    Args:
        a (int): Base
        n (int): Exponent value
        p (int): Divisor 

    Returns:
        int: a^n mod p
    """

    ret = 1
    while n > 0:
        if n % 2 == 1:
            ret = (ret * a) % p
        a = (a * a) % p
        n //= 2
    return ret

def legendre_num(a, p):
    """Gets the Legendre Number

    Args:
        a (int): Base
        p (int): Some prime

    Returns:
        int: Returns 1 if a is a quadratic residue modulo p 
    """

    return exp_mod(a, (p - 1) // 2, p)

def find_factor_base(N, B):
    """Generates factor base using Euler’s criterion to determine whether N is a quadratic residue mod p

    Args:
        N (int): Number to factor
        B (int): Upper bound to stop searching for primes

    Returns:
        List[int]: A list of primes in the factor base
    """

    if B < 2:
        return []
    
    prime_bool = [True for i in range(B + 1)]

    #start at first prime
    p = 2 

    factor_base = []

    #look for primes up until B
    while p < (B + 1): 
        if prime_bool[p]:
            #factor base only includes primes with legrendre num equal to 1
            if legendre_num(N,p) == 1: 
                    factor_base.append(p) 

            #mark off all multiples of p as being composite
            for j in range(2 * p, B + 1, p):
                prime_bool[j] = False
        p += 1
    return factor_base

def set_B(N, epsilon=0.1): 
    """Picks optimal bound B

    Args:
        N (int): Number to factor
        epsilon (float, optional): Some value close to 0. Defaults to 0.1

    Returns:
        int: Bound B for factor base
    """
    
    e = 2.71828183 
    #formula is optimal value for B
    #epsilon varies based on size of N
    B = e ** ((0.5 + epsilon) * sqrt(ceil(log(N)*log(log(N)))))
    return ceil(B)

def find_smooth_numbers_brute_force(factor_base, N, interval, extra_rows=1):
    """Finds B-smooth numbers by checking if the prime factors are in the factor base

    Args:
        factor_base (List[int]): Factor base
        N (int): Number to factor
        interval (int): Range of values to check for B-smooth numbers
        extra_rows (int, optional): Number of extra B-Smooth numbers to find. Defaults to 1

    Returns:
        List[int]: A list of B-smooth numbers
        List[int]: A list of x+n corresponding to each smooth number
        Dict[int:List[int]]: A dictionary mapping a smooth number to its prime factors
    """

    root = ceil(sqrt(N))
    x = [] 
    smooth_nums = []
    factors = {}


    #calculate (x + n)^ - N for n=0,1,2,3... in the interval
    for i in range(interval):
        
        f_i = int((root + i) * (root + i) - N)
        if f_i > N:
            f_i = f_i % N

        if f_i == 0:
            #this would create an infinite loop a few lines later
            continue
    
        f_i_original = f_i
        prime_factors = []
        for p in factor_base: 
            #p is a prime factor
            while f_i % p == 0:
                #divide by prime in factor base
                f_i //= p 
                prime_factors.append(p)
        
        if f_i == 1:
            smooth_nums.append(f_i_original)
            factors[f_i_original] = prime_factors
            x.append(i+root)
            #print every 10 smooth numbers
            if (len(smooth_nums) % 10 == 0):
                print('Smooth Number ' + str(len(smooth_nums))  + ':'+ str(f_i_original))
        #we have found enough smoother numbers
        if len(smooth_nums) >= (len(factor_base) + extra_rows): 
              break
    return smooth_nums, x, factors

def lower_register(orig, p):
    """Repeatedly divides a number by p 

    Args:
        orig (int): Some number
        p (int): Some prime

    Returns:
        int: The orginal number with factors of p divided
    """
    
    while (orig % p) == 0:
        orig //= p
    return orig

def lower_register_log(orig, log_p):
    return orig - log_p

def check_smooth(f, x, N, factor_base, smooth_nums, x_list, factors):
    """Checks if a number is smooth

    Args:
        f (int): Smooth number to check
        x (int): Square root of N
        N (int): Number to factor
        factor_base (List[int]): Factor base
        smooth_nums (List[int]): List of smooth numbers
        x_list (List[int]):  A list of x+n corresponding to each smooth number
        factors (Dict[int:List[int]]): A dictionary mapping a smooth number to its prime factors
        threshold (int): Number of extra smooth numbers
    """

    if f < 1:
        return
    if f == 1:
        #f(x) is smooth, where f(x) = (x)^2 - N
        x_list.append(x)
        f_orig = x * x - N
        smooth_nums.append(f_orig)

        #now find factors of f_orig
        prime_factors = []
        for p in factor_base:
            while f_orig % p == 0:
                f_orig //= p
                prime_factors.append(p)
        factors[x * x - N] = prime_factors

def check_smooth_log(f, x, N, factor_base, smooth_nums, x_list, factors, threshold):
    """Checks if a number is smooth, but log version

    Args:
        f (int): Smooth number to check
        x (int): Square root of N
        N (int): Number to factor
        factor_base (List[int]): Factor base
        smooth_nums (List[int]): List of smooth numbers
        x_list (List[int]):  A list of x+n corresponding to each smooth number
        factors (Dict[int:List[int]]): A dictionary mapping a smooth number to its prime factors
        threshold (int): Number of extra smooth numbers
    """

    if f <= 0:
        return
    elif f < threshold:

        #f(x) is smooth, where f(x) = (x)^2 - N
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

def find_smooth_numbers_TS(factor_base, N, extra_rows=1):
    """Finds B-smooth numbers with Tonelli Shanks algorithim

    Args:
        factor_base (int): Factor Base
        N (int): Number to factor
        extra_rows (int, optional): Number of extra B-Smooth numbers to find. Defaults to 1

    Returns:
        List[int]: A list of B-smooth numbers
        List[int]: A list of x+n corresponding to each smooth number
        Dict[int:List[int]]: A dictionary mapping a smooth number to its prime factors
    """

    root = ceil(sqrt(N))
    x_list = [] 
    smooth_nums = []
    factors = {}
    num_needed = len(factor_base) + extra_rows
    print("Number of smooth numbers needed is " + str(num_needed))

    solutions_mod_p = [tonelli_shanks(i, N) for i in factor_base]

    length = max(factor_base) * 10

    end=-1
    while len(smooth_nums) < num_needed:
        start=end+1
        end = start+length

        f_x = [(root + start + i) * (root + start + i) - N for i in range(length)]
       
        for i in range(len(factor_base)):
            p_i = factor_base[i]
            (x1, x2) = solutions_mod_p[i]
            #iterate through the f_x array dividing out p_i 

            #first val s.t. f_x[first_val_1] divisible by p_i
            first_val_1 = (0 - ((root + start) // (0 - p_i))  * p_i) - start - root + x1
            first_val_2 = (0 - ((root + start) // (0 - p_i))  * p_i) - start - root + x2
            if first_val_1 >= p_i:
                first_val_1 -= p_i
            if first_val_2 >= p_i:
                first_val_2 -= p_i
            
            for x_val in range(first_val_1, length, p_i):
                #divide out p
                f_x[x_val] = lower_register(f_x[x_val], p_i)
                
                #check that x_val+start is correct argument
                check_smooth(f_x[x_val], x_val+start+root, N, factor_base, smooth_nums, x_list, factors)
             

            for x_val in range(first_val_2, length, p_i):
                #divide out p
                f_x[x_val] = lower_register(f_x[x_val], p_i)
                check_smooth(f_x[x_val], x_val+start+root, N, factor_base, smooth_nums, x_list, factors)

    return smooth_nums, x_list, factors

def make_matrix(factor_base, smooth_nums, factors):
    """Builds matrix of exponents of prime factors of smooth numbers

    Args:
        factor_base (List[int]): Factor base
        smooth_nums (List[int]): List of smooth numbers
        factors (Dict[int:List[int]]): A dictionary mapping a smooth number to its prime factors

    Returns:
        List[List[int]]: A 2d matrix mod 2
        List[List[int]]: A 2d matrix 
    """
    
    #binary version
    matrix = [] 
    #non-binary version
    matrix_nb = [] 

    for n in smooth_nums:
        exp_vector = [0 for i in range(len(factor_base))]
        exp_vector2 = [0 for i in range(len(factor_base))]
        prime_factors = factors[n]

        for i in range(len(factor_base)):
            if factor_base[i] in prime_factors:
                #get exponenet vector of prime factors mod 2
                exp_vector[i] = (exp_vector[i] + prime_factors.count(factor_base[i])) % 2 
                exp_vector2[i] = prime_factors.count(factor_base[i])
        matrix.append(exp_vector)
        matrix_nb.append(exp_vector2)
    return matrix, matrix_nb

def matrix_transpose(matrix):
    """Tranposes a matrix

    Args:
        matrix (List[List[int]]): Some matrix

    Returns:
        List[List[int]]: The tranposed versioned of the matrix
    """

    original_num_rows = len(matrix)
    original_num_cols = len(matrix[0])
    new_matrix = [[0 for i in range(original_num_rows)] for j in range(original_num_cols)]
    for i in range(original_num_cols):
        for j in range(original_num_rows):
            new_matrix[i][j] = matrix[j][i]
    return new_matrix

def find_left_null_vectors(matrix):
    """Finds the left null space vector of a matrix

    Args:
        matrix (List[List[int]]): Some matrix

    Returns:
        List[int]: A vector representing the left null space
    """

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

def find_squared_congruence(null_vector, x, factor_base, matrix_nb, N):
    """Caulates a and b for a^2 = b^2 mod n from left null space

    Args:
        null_vector (List[int]): Vector representing left null space
        x (List[int]): A list of x+n corresponding to each smooth number 
        factor_base (List[int]): Factor Base
        matrix_nb (List[List[int]]): Matrix representing prime factors for each smooth number
        N (int): Number to factor

    Returns:
        (int, int): a and b
    """

    combo_rows = []
    for i in range(len(null_vector)):
        if null_vector[i] == 1:
            combo_rows.append(i)
    x_nums = [x[i] for i in combo_rows]
    
    a = 1
    b = 1
    #this is an exponent vector
    smooth_vector = [0 for i in factor_base] 

    for i in combo_rows:
        for j in range(len(factor_base)):
            #gets the exponent power for each prime
            smooth_vector[j] += matrix_nb[i][j] 
  
    for i in range(len(smooth_vector)):
        f = 1
        for j in range(smooth_vector[i] >> 1):
            f = f * (factor_base[i])
        a *= f 

    for n in x_nums:
        b *= n
    
    if ((a * a) % N ) != ((b * b) % N):
        print("Error, problem!")
    return (a, b)   

def check_prime(n, iter=10):
    """Runs Miller Rabin to check if a number is prime

    Args:
        n (int): Number to check if prime
        iter (int, optional): Number of iterations. Defaults to 10

    Returns:
        bool: True if prime and False otherwise
    """

    if (exp_mod(2, n-1, n) != 1):
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

def tonelli_shanks(p, n):
    """Tonellis Shanks algorithim to find solutions for x^2 = n mod p

    Args:
        p (int): Some prime
        n (int): Some number

    Returns:
        (int, int): 2 solutions x and p - x
    """

    if legendre_num(n, p) != 1:
        print("not a square (mod p)")

    #now going to find s and e s.t. s * 2^e = p-1
    s = p - 1
    e = 0
    while (s % 2) == 0:
        e += 1
        s //= 2
    if e == 1:
        x_1 = exp_mod(n, (p + 1) // 4, p)
        x_2 = (p - x_1) % p
        return (x_1, x_2)

    #now finding q s.t. q^{(p-1)/2} is congruent to -1 mod p
    q = 2
    for q_test in range(2, p, 1):
        if legendre_num(q_test, p) == p-1:
            q = q_test
            break

    #initializing some variables
    x = exp_mod(n, (s+1)//2, p)
    b = exp_mod(n, s, p)              
    g = exp_mod(q, s, p)
    r = e  

    #now iterating until we find a solution
    while (b-1) % p != 0:
        b_sq = (b * b) % p
        i = 1
        while (i < r):
            if (b_sq - 1) % p == 0:
                break
            b_sq = (b_sq * b_sq) % p
            if (i + 1 < r):
                i+= 1
            else:
                break
        
        a = exp_mod(g, 1 << (r - i - 1), p)
        x = (x * a) % p
        g = (a * a) % p
        b = (b * g) % p
        r = i
    return (x, (p - x) % p)

def factor(N, epsilon):
    """Full quadratic sieve algorithim that factors a number

    Args:
        N (int): Number to factor
        epsilon (float): Epsilon to pick bound B

    Returns:
        List[int]: The times to factor N in milliseconds
    """

    start_time_ms = int(time.time() * 1000)

    if check_prime(N):
        print("N is very likely prime")
        return

    B = set_B(N, epsilon)

    print('Bound B: ' + str(B))
    
    factor_base = find_factor_base(N, B)

    time_0 = int(time.time() * 1000)
    print('Finding smooth numbers ...')
    extra_rows = int(20 + log(N))

    #use brute force if N is small otherwise use Tonelli Shanks method
    if (N < ( 10 ** 15 )):
        smooth_nums, x, factors = find_smooth_numbers_brute_force(factor_base, N, B**3, extra_rows)
    else:
        smooth_nums, x, factors = find_smooth_numbers_TS(factor_base, N, extra_rows)
        if smooth_nums == None:
            smooth_nums, x, factors = find_smooth_numbers_brute_force(factor_base, N, B**3, extra_rows)
    found_smooth_time_ms = (time.time() * 1000) - time_0

    if len(smooth_nums) < len(factor_base):
        print('Not enough B-smooth numbers were found')
        return
    else:
        print('Found smooth numbers in ' + str(int(found_smooth_time_ms)) + "ms")

    #below, matrix is all mod 2 entries, while matrix_nb has any positive int for each entry
    matrix, matrix_nb = make_matrix(factor_base, smooth_nums, factors) 

    time_1 = int(time.time() * 1000)
    print('Running Gaussian elimination ... ')
    
    null_combos = find_left_null_vectors(matrix)
    elim_time = int(int(time.time() * 1000) - time_1)
    print("Finished Gaussian elimination in: " + str(int(elim_time)) + "ms")

    num_solutions = len(null_combos)

    for s in range(num_solutions):
        [a, b] = find_squared_congruence(null_combos[s], x, factor_base, matrix_nb, N)
        possible_factor = abs(gcd(a-b, N))
        if possible_factor == 1 or possible_factor == N:
            continue
        else:
            print("N = " + str(possible_factor) + " * " + str(int(N//possible_factor)))
            total_time = int(time.time() * 1000) - start_time_ms
            print("Total time: " + str(int(total_time)) + "ms")
            return [possible_factor, int(N//possible_factor)]
        
    #if the code reaches this point, no nontrivial factors have been found
    #retry with higher epsilon
    print("Could not find a nontrivial factor, retrying to factor")
    return factor(N, epsilon * 10)

def main():
    N = 16921456439215439701
    return factor(N, epsilon=0.1)

if __name__ == "__main__":
    print(main())

