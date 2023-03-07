"""
Output:
Solution:
[7.0, 5.000000000000001, 4.000000000000001, 6.0, 4.000000000000001]
Trigonometric matrix:
[[-11, -9, 0, 0, 0], [5, -15, -2, 0, 0], [0, -8, 11, -3, 0], [0, 0, 6, -15, 4], [0, 0, 0, 3, 6]]
[-122, -48, -14, -50, 42]
"""


def solve(m, b):
    y = [None] * len(m)
    alpha = [None] * len(m)
    beta = [None] * len(m)
    for i in range(len(m)):
        if i == 0:
            y[i] = m[i][0]
            alpha[i] = -1 * m[i][1] / y[i]
            beta[i] = b[i] / y[i]
        elif i == len(m) - 1:
            y[i] = m[i][1] + m[i][0] * alpha[i-1]
            beta[i] = (b[i] - m[i][0] * beta[i-1]) / y[i]
        else:
            y[i] = m[i][1] + m[i][0] * alpha[i-1]
            alpha[i] = -1 * m[i][2] / y[i]
            beta[i] = (b[i] - m[i][0] * beta[i-1]) / y[i]
    x = [0] * len(m)
    for i in range(len(m)):
        if i == 0:
            x[len(m)-i-1] = beta[len(m)-i-1]
        else:
            x[len(m)-i-1] = alpha[len(m)-i-1] * x[len(m)-i] + beta[len(m)-i-1]
    return x


def trig_matrix(A):
    m = []
    for i in range(len(A)):
        m.append([0]*len(A))
    for i in range(len(A)):
        if i == 0:
            m[i][0] = A[i][0]
            m[i][1] = A[i][1]
        elif i == len(A)-1:
            m[i][i-1] = A[i][0]
            m[i][i] = A[i][1]
        else:
            m[i][i-1] = A[i][0]
            m[i][i] = A[i][1]
            m[i][i+1] = A[i][2]
    return m


A = [[-11,-9],
     [5,-15,-2],
     [-8,11,-3],
     [6,-15,4],
     [3,6]]

b = [-122, -48, -14, -50, 42]

x = solve(A, b)
print("Solution:")
print(x)
print("Trigonometric matrix:")
print(trig_matrix(A))
print(b)
