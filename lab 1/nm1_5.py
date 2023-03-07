"""
Output:
iter: 30
[[-13.501038338364921, -0.5056959789742247, -4.352049331191065], [-6.732306350240983e-05, 3.9255012743783517, -4.603892322209252], [-1.7498775222900943e-05, 2.163984641208537, 6.575537063986573]]
roots:
[-13.501038338364921, [5.250519169182462, 2.864800141979823], [5.250519169182462, -2.864800141979823]]
"""

import math
import copy


def solve(A):
    iter = 0
    Q,R = find_QR(A)
    A_ = prois(Q, R)
    while error(A, A_, err) != True:
        Q,R = find_QR(A)
        A = A_
        A_ = prois(Q, R)
        iter += 1
    print('iter:',iter)
    return A_


def solve_roots(A):
    res = [A[0][0]]
    x = A[1][1]
    y = A[2][2]
    z = A[1][2] * A[2][1]
    D = (x+y)*(x+y)-4*(x*y-z)
    res.append([(x+y)/2, math.sqrt(-D)/2])
    res.append([(x+y)/2, -math.sqrt(-D)/2])
    print(res)


def display(m, t):
    print(t)
    n = 5
    template = '{:.' + str(n) + 'f}'
    for row in m:
        for num in row:
            print(str(template.format(num)), end='  ')
        print()
    print('')


def m_zero(n):
    res = []
    for i in range(n):
        r = []
        for j in range(n):
            r.append(0)
        res.append(r)
    return res


def prois(U, L):
    n = len(U)
    R = m_zero(n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                R[i][j] += L[i][k] * U[k][j]
    return R


def find_housholder(A, it):
    v = [0] * len(A)
    for i in range(it, len(A)):
        v[i] = A[i][it]
    s = 0
    for i in range(len(A)):
        s += math.pow(v[i],2)
    if A[it][it] < 0:
        v[it] -= math.sqrt(s)
    elif A[it][it] > 0:
        v[it] += math.sqrt(s)
    else:
        print('eq 0 house')
    H = copy.deepcopy(A)
    dim = 0
    for i in range(len(v)):
        dim += math.pow(v[i],2)
    for i in range(len(A)):
        for j in range(len(A)):
            if i == j:
                H[i][j] = 1 - 2 * v[i] * v[j] / dim
            else:
                H[i][j] = 0 - 2 * v[i] * v[j] / dim
    return H


def find_QR(A):
    R_ = copy.deepcopy(A)
    Q_ = None
    for i in range(len(R_) - 1):
        H = find_housholder(R_, i)
        if Q_ == None:
            Q_ = H
        else:
            Q_ = prois(H, Q_)
        R_ = prois(R_, H)
    return Q_, R_


def error(A, A_, err):
    s = 0
    for i in range(len(A)):
        for j in range(len(A)):
            if j == 0 and i > j:
                s += math.pow(A_[i][j], 2)
    if math.sqrt(s) < err:
        return True
    return False


err = 10**-4

A = [[3,-7,-1],
    [-9,-8, 7],
    [ 5, 2, 2]]

matrix = solve(A)
print(matrix)
print("roots:")
solve_roots(matrix)
