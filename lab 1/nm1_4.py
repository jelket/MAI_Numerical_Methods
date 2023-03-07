"""
Output:
Matrix:
[[2.073768569795485, 0.008523674778537393, -0.000995859209226302], [0.00852367477857039, -19.363120113404072, -0.0014381164706058769], [-0.0009958592092440483, -0.001438116470625861, -3.7106484563913593]]
Matrix eigenvalues:
[-19.363120113404072, -3.7106484563913593, 2.073768569795485]
"""

import copy
import math


def jacobi(A, err):
    A_ = copy.deepcopy(A)
    num_of_it = 0
    while True:
        i, j = find_max(A)
        P = math.pi / 4
        if A[i][i] - A[j][j] != 0:
            P = 2 * A[i][j] / (A[i][i] - A[j][j])
        c = math.cos(math.atan(P) / 2)
        s = math.sin(math.atan(P) / 2)
        rotate = rotate_matrix(A,s,c,i,j)
        A_ = prois(transpose(rotate),prois(A,rotate))
        num_of_it += 1
        A = copy.deepcopy(A_)
        if error(A_, err):
            break
    return A, num_of_it


def find_max(A):
    m, ib, jb = None, None, None
    for i in range(len(A)):
        for j in range(len(A)):
            if i < j:
                if m == None or abs(A[i][j]) > m:
                    m = A[i][j]
                    ib, jb = i, j
    return ib, jb


def error(A_, err):
    s = sum([sum([math.pow(A_[i][j], 2) for j in range(i+1,len(A_))]) for i in range(len(A_))])
    return False if math.sqrt(s) > err else True


def transpose(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A))]


def rotate_matrix(A,s,c,i,j):
    res = copy.deepcopy(A)
    for k in range(len(A)):
        for l in range(len(A)):
            if k == l:
                res[k][l] = 1
            else:
                res[k][l] = 0
    res[i][i] = c
    res[i][j] = -s
    res[j][i] = s
    res[j][j] = c
    return res


def prois(U, L):
    n = len(U)
    return [[sum([L[i][k]*U[k][j] for k in range(n)]) for j in range(n)] for i in range(n)]


def get_self_number(A):
    return sorted([A[i][i] for i in range(len(A))])


err = 0.01
A = [[-7, 4, 5], [4, -6, -9], [5, -9, -8]]
matrix, cnt = jacobi(A, err)
print("Matrix:")
print(matrix)
print("Matrix eigenvalues:")
print(get_self_number(matrix))
