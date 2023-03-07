import numpy as np
from prettytable import PrettyTable


def runge_kutta(f, xa, xb, ya, y1a, h):
    n = int((xb - xa) / h)
    x = xa
    y = ya
    z = y1a
    x_res = [x]
    y_res = [y]
    z_res = [z]
    for i in range(1, n + 1):
        k1 = h * z
        l1 = h * f(x, y, z)
        k2 = h * (z + 0.5 * l1)
        l2 = h * f(x + 0.5 * h, y + 0.5 * k1, z + 0.5 * l1)
        k3 = h * (z + 0.5 * l2)
        l3 = h * f(x + 0.5 * h, y + 0.5 * k2, z + 0.5 * l2)
        k4 = h * (z + l3)
        l4 = h * f(x + h, y + k3, z + l3)
        x = xa + i * h
        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        z += (l1 + 2 * l2 + 2 * l3 + l4) / 6
        x_res.append(x)
        y_res.append(y)
        z_res.append(z)
    return (x_res, y_res), z_res


def finite_difference_method(a1, b1, h, alpha_0, alpha_1, beta_0, beta_1, A, B, fi=None):
    x = [a1]
    a = []
    b = []
    c = []
    d = []
    n = round((b1 - a1) / h)
    a.append(0)
    b.append(-2 / (h * (2 - p(a1) * h)) + q(a1) * h /
             (2 - p(a1) * h) + alpha_0 / beta_0)
    c.append(2 / (h * (2 - p(a1) * h)))
    d.append(A / beta_0 + h * f(a1) / (2 - p(a1) * h))
    x.append(x[0] + h)
    if fi:
        fi.write(f'{b[0]}y_0 + {c[0]}y_1 = {d[0]}\n')

    for i in range(1, n):
        a.append(1 / h ** 2 - p(x[i]) / (2 * h))
        b.append(-2 / h ** 2 + q(x[i]))
        c.append(1 / h ** 2 + p(x[i]) / (2 * h))
        d.append(f(x[i]))
        x.append(x[i] + h)
        if fi:
            fi.write(f'{a[i]}y_{i - 1} + {b[i]}y_{i} + {c[i]}y_{i+1} = {d[i]}\n')
    a.append(-2 / (h * (2 + p(x[n]) * h)))
    b.append(2 / (h * (2 + p(x[n]) * h)) - q(x[n]) * h /
             (2 + p(x[n]) * h) + alpha_1 / beta_1)
    c.append(0)
    d.append(B / beta_1 - h * f(x[n]) / (2 + p(x[n]) * h))
    if fi:
        fi.write(f'{a[-1]}y_{len(a) - 2} + {b[-1]}y_{len(a) - 1} = {d[-1]}\n')
    y = tma(a, b, c, d, len(a))
    return x, y


def func(x, y, y1):
    return (x * y - 2 * y1) / x


def orig_func(x):
    return np.exp(x) / x


def p(x):
    return 2 / x


def q(x):
    return -1


def f(x):
    return 0


def tma(a, b, c, d, shape):
    p = [-c[0] / b[0]]
    q = [d[0] / b[0]]
    x = [0] * (shape + 1)
    for i in range(1, shape):
        p.append(-c[i] / (b[i] + a[i] * p[i - 1]))
        q.append((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]))
    for i in reversed(range(shape)):
        x[i] = p[i] * x[i + 1] + q[i]
    return x[:-1]


def der_one(xi, yi, x):
    i = 0
    while xi[i + 1] < x - 1e-7:
        i += 1
    return (yi[i + 1] - yi[i]) / (xi[i + 1] - xi[i])


def next_n(cur_n, prev_n, ans_cur, ans_prev, alpha1, beta1, B, b, fi=None):
    num1 = beta1 * der_one(ans_cur[0], ans_cur[1], b)
    num2 = beta1 * der_one(ans_prev[0], ans_prev[1], b)
    num3 = alpha1 * ans_prev[1][len(ans_prev[0]) - 1]
    num4 = alpha1 * ans_cur[1][len(ans_cur[0]) - 1] + num1 - B
    num5 = alpha1 * ans_cur[1][len(ans_cur[0]) - 1] + num1 - num3 - num2
    if fi:
        fi.write(f'{round(cur_n, 5)} - ({round(cur_n, 5)} - {round(prev_n, 5)}) / '
                 f'({round(alpha1 * ans_cur[1][len(ans_cur[0]) - 1] + num1, 5)} - {round(num3 - num2, 5)}) * '
                 f'({round(alpha1 * ans_cur[1][len(ans_cur[0]) - 1] + num1, 5)} - {round(B, 5)})'
                 f' = {round(cur_n - num4 * (cur_n - prev_n) / num5, 5)}\n')
    return cur_n - num4 * (cur_n - prev_n) / num5


def shooting_method(a, b, h, eps, f, alpha0, alpha1, beta0, beta1, A, B, fi=None):
    table = PrettyTable(['j', 'n_j', 'y', '|Phi(n_j)|'])

    n_prev = 1
    n_cur = 0.8
    if fi:
        fi.write(f'|{alpha1} * y({a}, n_k, {(A - alpha0 * n_prev) / beta0}) '
                 f'+ {beta1} * y\'({a}, n_k, {(A - alpha0 * n_prev) / beta0})) - {B}| <= eps\n')
        fi.write(f'eps = {eps}\n')
    ans_prev = runge_kutta(f, a, b, n_prev, (A - alpha0 * n_prev) / beta0, h)[0]
    ans_cur = runge_kutta(f, a, b, n_cur, (A - alpha0 * n_cur) / beta0, h)[0]
    table.add_row([0, n_prev, ans_prev[1][-1], abs(alpha1 * ans_prev[1][len(ans_prev[0]) - 1] +
                                                   beta1 * der_one(ans_prev[0], ans_prev[1], b) - B)])
    table.add_row([1, n_cur, ans_cur[1][-1], abs(alpha1 * ans_cur[1][len(ans_cur[0]) - 1] +
                                                 beta1 * der_one(ans_cur[0], ans_cur[1], b) - B)])
    i = 2
    while abs(alpha1 * ans_cur[1][len(ans_cur[0]) - 1] +
              beta1 * der_one(ans_cur[0], ans_cur[1], b) - B) > eps:
        if i == 2 or i == 3:
            n = next_n(n_cur, n_prev, ans_cur, ans_prev, alpha1, beta1, B, b, fi)
        else:
            n = next_n(n_cur, n_prev, ans_cur, ans_prev, alpha1, beta1, B, b)
        n_prev = n_cur
        n_cur = n
        ans_prev = ans_cur
        ans_cur = runge_kutta(f, a, b, n_cur, (A - alpha0 * n_cur) / beta0, h)[0]
        table.add_row([i, n_cur, ans_cur[1][-1], abs(alpha1 * ans_cur[1][len(ans_cur[0]) - 1] + \
                                                     beta1 * der_one(ans_cur[0], ans_cur[1], b) - B)])
        i += 1
    if fi:
        fi.write(f'\n{str(table)}\n')
    return ans_cur


def print_result_table(f, name, res, runge_y, p):
    f.write(f'\n{name}\n')
    table = PrettyTable(['x', 'y', 'y_src', 'eps_k',
                         'Romberg\'s method error'])
    for x, y, yr in zip(*res, runge_y):
        tmp = orig_func(x)
        table.add_row([round(x, 2), round(y, 5), round(tmp, 5), abs(round(y - tmp, 5)),
                       abs(round(y + (y - yr) / (0.5 ** p - 1) - tmp, 8))])
    f.write(f'\n{str(table)}\n')


def main():

    a = 1
    b = 2
    alpha0 = 0
    alpha1 = 1.5
    beta0 = 1
    beta1 = 1
    y0 = 0
    y10 = np.e ** 2
    step = 0.1
    eps = 0.001
    f = open("Output", 'w')

    res1 = shooting_method(a, b, step, eps, func, alpha0, alpha1,
                           beta0, beta1, y0, y10, f)
    res2 = shooting_method(a, b, step / 2, eps, func, alpha0, alpha1,
                           beta0, beta1, y0, y10)
    h_half = [y for x, y in zip(res2[0], res2[1]) if x in res1[0]]
    print_result_table(f, 'Shooting method', res1, h_half, 2)

    res1 = finite_difference_method(a, b, step, alpha0, alpha1,
                                    beta0, beta1, y0, y10, f)
    res2 = finite_difference_method(a, b, step / 2, alpha0, alpha1,
                                    beta0, beta1, y0, y10)
    h_half = [y for x, y in zip(res2[0], res2[1]) if x in res1[0]]
    print_result_table(f, 'Finite difference method', res1, h_half, 2)

    f.close()


if __name__ == "__main__":
    main()
