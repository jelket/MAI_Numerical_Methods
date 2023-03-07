"""
Output:
1.0c2 + 4.0c3 + 1.0c4 = -0.6962

0.0
0.5
0.86603
1.0
0.86603
x = [0.0, 1.0, 2.0, 3.0, 4.0]
y = [0.0, 0.5, 0.86603, 1.0, 0.86603]

a = [0.0, 0.5, 0.86603, 1.0]
b = [0.52409375, 0.45181249999999995, 0.26674624999999996, -0.01879749999999998]
c = [0, -0.07228125000000005, -0.11278499999999991, -0.17275875000000007]
d = [-0.024093750000000014, -0.013501249999999956, -0.019991250000000054, 0.057586250000000026]

f(1.5) = 0.7061
"""

import numpy as np
import matplotlib.pyplot as plt


def tma(matrix, d, shape):
    a, b, c = zip(*matrix)
    p = [-c[0] / b[0]]
    q = [d[0] / b[0]]
    x = [0] * (shape + 1)
    for i in range(1, shape):
        p.append(-c[i] / (b[i] + a[i] * p[i - 1]))
        q.append((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]))
    for i in reversed(range(shape)):
        x[i] = p[i] * x[i + 1] + q[i]
    return x[:-1]


def spline(x, y):
    size = len(x)
    h = [x[i] - x[i - 1] for i in range(1, size)]
    mtrx = [[0, 2 * (h[0] + h[1]), h[1]]]
    b = [3 * ((y[2] - y[1]) / h[1] - (y[1] - y[0]) / h[0])]
    for i in range(1, size - 3):
        tmp = [h[i], 2 * (h[i] + h[i + 1]), h[i + 1]]
        mtrx.append(tmp)
        b.append(3 * ((y[i + 2] - y[i + 1]) / h[i + 1] - (y[i + 1] - y[i]) / h[i]))
    mtrx.append([h[-2], 2 * (h[-2] + h[-1]), 0])
    b.append(3 * ((y[-1] - y[-2]) / h[-1] - (y[-2] - y[-3]) / h[-2]))
    print(f'{mtrx[i][0]}c{i+1} + {mtrx[i][1]}c{i+2} + {mtrx[i][2]}c{i+3} = {round(b[i],4)}\n')
    c = tma(mtrx, b, size - 2)
    a = []
    b = []
    d = []
    c.insert(0, 0)
    for i in range(1, size):
        a.append(y[i - 1])
        if i < size - 1:
            d.append((c[i] - c[i - 1]) / (3 * h[i - 1]))
            b.append((y[i] - y[i - 1]) / h[i - 1] -
                     h[i - 1] * (c[i] + 2 * c[i - 1]) / 3)
    b.append((y[-1] - y[-2]) / h[-1] - 2 * h[-1] * c[-1] / 3)
    d.append(-c[-1] / (3 * h[-1]))
    return a, b, c, d


def polyval(x0, x, k, coef):
    a, b, c, d = coef
    tmp = (x0 - x[k])
    return a[k] + b[k] * tmp + c[k] * tmp**2 + d[k] * tmp**3


def pol(x, x_test, coef):
    k = 0
    for i, j in zip(x, x[1:]):
        if i <= x_test <= j:
            break
        k += 1
    return polyval(x_test, x, k, coef)


def main():

    x = [0.0, 1.0, 2.0, 3.0, 4.0]
    y = [0.0, 0.5, 0.86603, 1.0, 0.86603]
    x_test= 1.5

    coef = spline(x, y)

    x1 = np.linspace(x[0], x[-1], 50)
    y1 = [pol(x, i, coef) for i in x1]

    plt.plot(x1, y1, color='b')
    plt.scatter(x, y, color='r')
    plt.show()

    res = pol(x, x_test, coef)

    for i in x:
        print(pol(x, i, coef))

    print(f'x = {x}\ny = {y}\n\n')
    print(f'a = {coef[0]}\nb = {coef[1]}\nc = {coef[2]}\nd = {coef[3]}\n')
    print(f'\nf({x_test}) = {round(res, 4)}\n')


if __name__ == "__main__":
    main()
