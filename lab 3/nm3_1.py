"""
Output:
a-Newton) 0.00016011185300468522
a-Lagrange) 0.00016011185300468522
b-Newton) 0.0002610156283867626
b-Lagrange) 0.00026101562838687364
"""


from math import pi, tan, sin, cos


def f(x):
    return sin(x)


x_a = [0.1*pi, 0.2*pi, 0.3*pi, 0.4*pi]
y_a = [f(_) for _ in x_a]

x_b = [0.1*pi, pi/6, 0.3*pi, 0.4*pi]
y_b = [f(_) for _ in x_b]

x_star = pi/4
y_star = f(x_star)


def divided_diff(x, y):
    n = len(y)
    coef = []
    for i in range(len(x)):
        r = [y[i]]
        for j in range(len(x)-1):
            r.append(0)
        coef.append(r)
    for j in range(1,n):
        for i in range(n-j):
            coef[i][j] = (coef[i+1][j-1] - coef[i][j-1]) / (x[i+j]-x[i])
    return coef


def newton(X, x, y):
    coef = divided_diff(x, y)
    res, cof = coef[0][0], []
    for i in range(1,len(coef)):
        cof.append(coef[0][i])
    for i in range(len(cof)):
        for j in range(i+1):
            cof[i] *= (X - x[j])
        res += cof[i]
    return res


def lagrange(X,x,y):
    res = 0
    for i in range(len(x)):
        f_i = y[i]
        for j in range(len(y)):
            if j != i:
                f_i *= (X - x[j]) / (x[i] - x[j])
        res += f_i
    return res


print('a-Newton)', abs(f(x_star) - newton(x_star, x_a,y_a)))
print('a-Lagrange)', abs(f(x_star) - lagrange(x_star, x_a,y_a)))

print('b-Newton)', abs(f(x_star) - newton(x_star, x_b,y_b)))
print('b-Lagrange)', abs(f(x_star) - lagrange(x_star, x_b,y_b)))
