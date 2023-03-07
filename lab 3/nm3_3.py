"""
Output:
n = 2 error = 0.12320284825395433
n = 3 error = 0.018643313866991175
"""

from math import pi, tan, sin, pow, sqrt
import copy


def make_spline_matrix(x, y, n, err=10**-5, text=False):
    A, b, n = [], [], n+1
    for i in range(n):
        r = []
        for j in range(n):
            if i == 0 and j == 0:
                r.append(len(x))
            else:
                r.append(sum(map(lambda a: pow(a,i+j),x)))
        A.append(r)
        b.append(sum(map(lambda a,b: pow(a,i) * b,x,y)))
    if text:
        print('Matrix A and b')
        for i in range(len(A)):
            for j in range(len(A[i])):
                print('{:.2f}'.format(A[i][j]), end='  ')
            print('['+'{:.2f}'.format(b[i])+']')
    a_,a = [None] * len(A),[0] * len(A)
    while True:
        for i in range(len(A)):
            s = 0
            for j in range(len(A)):
                if j < i:
                    s += A[i][j] * a_[j]
                elif i != j:
                    s += A[i][j] * a[j]
            a_[i] = (b[i] - s) / A[i][i]
        if sqrt(sum(map(lambda a,b: pow(a - b,2),a,a_))) < err:
            break
        a = copy.copy(a_)
    if text:
        print('f(x) = ', float('{:.4f}'.format(a_[0])), end='')
        for i in range(1, len(a_)):
            print(' +', str(float('{:.4f}'.format(a_[i]))) + '*x^' + str(i),end='')
        print()

    resx,resy = [],[]
    start = x[0]
    while start < x[-1]:
        resx.append(start)
        resy.append(sum([a_[j] * pow(start, j) for j in range(len(a_))]))
        start += 0.1
    yy = [sum([a_[j] * pow(num, j) for j in range(len(a_))]) for num in x]
    print('n =', n-1, 'error =', sqrt(sum(list(map(lambda a, b: pow(a - b,2), yy, y)))))
    return resx, resy


x_star = 2.666667

x = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
y = [-0.5, 0.0, 0.5, 0.86603, 1.0, 0.86603]


x_2, y_2 = make_spline_matrix(x, y, 2)
x_3, y_3 = make_spline_matrix(x, y, 3)
