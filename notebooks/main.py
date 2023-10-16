import sympy as sp


# define parameter
t = sp.Symbol('t')

my_polynomial = t**2 - 1

roots = sp.roots(my_polynomial)

p = sp.Poly(my_polynomial, t, domain='RR')
coeffs = p.all_coeffs()
deg = p.degree()

# companion matrix
companion_matrix = sp.zeros(deg)
for i in range(deg - 1):
    companion_matrix[i + 1, i] = 1
for i in range(deg):
    companion_matrix[i, -1] = -1 * coeffs[::-1][i]

newton_sums = []
for k in range(2*deg - 1):
    V_k = companion_matrix ** k
    newton_sums.append(V_k.trace())


# Get Newton sums
# https://brilliant.org/wiki/newtons-identities/
#n = sp.symbols('n')
#newton_sums2 = [deg, -1*coeffs[1], coeffs[1]**2 - 2*coeffs[2], -1*coeffs[1]**3 + 3*coeffs[1]*coeffs[2] - 3*coeffs[3]]
#print(newton_sums)

# Construct the Hankel matrix
hankel_matrix = sp.zeros(deg)
for i in range(deg):
    for j in range(deg):
        hankel_matrix[i, j] = newton_sums[i+j]

eigenvalues = hankel_matrix.eigenvals()

# Count the positive and negative eigenvalues
positive_eigenvalues = sum(1 for eigenvalue in eigenvalues if eigenvalue > 0)
negative_eigenvalues = sum(1 for eigenvalue in eigenvalues if eigenvalue < 0)

# Calculate the signature
signature = positive_eigenvalues - negative_eigenvalues
rank = hankel_matrix.rank()

print(f"Rank: {rank}")
print(f"Sign: {signature}")
print('---------------------')
print('Roots: ')
print(roots)
