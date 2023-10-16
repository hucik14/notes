import sympy as sp


# define parameter
t = sp.Symbol('t')

my_polynomial = t**2 + 1

roots = sp.roots(my_polynomial)

p = sp.Poly(my_polynomial, t, domain='RR')
coeffs = p.all_coeffs()
deg = p.degree()

# Get Newton sums
# https://brilliant.org/wiki/newtons-identities/
n = sp.symbols('n')
newton_sums = [deg, -1*coeffs[1], coeffs[1]**2 - 2*coeffs[2]]


# Construct the Hankel matrix
hankel_matrix = sp.zeros(deg)
for i in range(deg):
    for j in range(deg):
        hankel_matrix[i, j] = newton_sums[i + j]

print(hankel_matrix)
