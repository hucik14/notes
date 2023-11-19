import sympy as sp
import numpy as np
import scipy


def symmetric_gaussian_elimination(A):
    n = A.shape[0]
    S = np.eye(n)  # Initialize S as the identity matrix

    for i in range(n - 1):
        for j in range(i + 1, n):
            # Compute Givens rotation parameters
            theta = np.arctan2(A[j, i], A[i, i])
            c = np.cos(theta)
            s = np.sin(theta)

            # Construct the Givens rotation matrix
            G = np.eye(n)
            G[i, i] = c
            G[j, j] = c
            G[i, j] = -s
            G[j, i] = s

            # Apply the Givens rotation to A and S
            A = G.T @ A @ G
            S = S @ G

    return S


# define parameter
t = sp.Symbol('t')

my_polynomial = t**2 - 1 - 2*t**3 + 3*t**2 - 4*t + 5

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
hermite_matrix = sp.zeros(deg)
for i in range(deg):
    for j in range(deg):
        hermite_matrix[i, j] = newton_sums[i+j]

# take hermite matrix as numpy matrix
hermite_np = np.zeros((deg, deg))
for i in range(deg):
    for j in range(deg):
        hermite_np[i, j] = newton_sums[i+j]

# Perform Cholesky decomposition
#L = np.linalg.cholesky(hermite_np)
#print(L)

# Perform LU decomposition with pivoting and unit lower triangular L
#P, L, U = scipy.linalg.lu(hermite_np, permute_l=True)

# Construct the diagonal matrix D from U
S = symmetric_gaussian_elimination(hermite_np)
D = S.transpose() @ hermite_np @ S
print(hermite_np)
print(D)

diag = np.diag(D)

# Count the positive and negative elements of D
positive_elements = sum(1 for val in diag if val > 0)
negative_elements = sum(1 for val in diag if val < 0)
sign = positive_elements - negative_elements


eigenvalues_dict = hermite_matrix.eigenvals()
eigenvalues = [key for key, value in eigenvalues_dict.items() for _ in range(value)]

# Count the positive and negative eigenvalues
positive_eigenvalues = sum(1 for eigenvalue in eigenvalues if eigenvalue > 0)
negative_eigenvalues = sum(1 for eigenvalue in eigenvalues if eigenvalue < 0)

# Calculate the signature
signature = positive_eigenvalues - negative_eigenvalues
rank = hermite_matrix.rank()

print('------------')
print('------------')
print(f"Rank: {rank}")
print(f"Sign: {sign}")
print(f"Sign: {signature} (using eigenvals)")
print(f'Roots: {roots}')
print('------------')
