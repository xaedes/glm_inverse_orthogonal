import sympy as sp

a,b,c,d,e,f = sp.symbols("a b c d e f")
orthogonal_matrix = sp.ImmutableDenseMatrix([
    [a,0,0,d],
    [0,b,0,e],
    [0,0,c,f],
    [0,0,0,1]
])

print(orthogonal_matrix)
print(orthogonal_matrix.inverse())

# Output:
# Matrix([[a, 0, 0, d], [0, b, 0, e], [0, 0, c, 0], [0, 0, 0, 1]])
# Matrix([[1/a, 0, 0, -d/a], [0, 1/b, 0, -e/b], [0, 0, 1/c, 0], [0, 0, 0, 1]])
