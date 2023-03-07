import numpy as np
from prettytable import PrettyTable


def foo(x, y, y1):
    return np.sin(3 * x) - y


def orig_foo(x):
    return np.cos(x) + 11 / 8 * np.sin(x) - np.sin(3 * x) / 8


def euler(f, xa, xb, ya, y1a, h, fl=None):
    if fl:
        fl.write(f'Euler:\n\n')
        fl.write(f'z = y\' \n'
                 f'z\' = sin(3x) - y \n'
                 f'[{xa}, {xb}]\n'
                 f'x0 = {xa}, y0 = {ya}, z0 = {y1a}\n\n')
    n = int((xb - xa) / h)
    x = xa
    y = ya
    x_res = [x]
    y_res = [y]
    deltaYk = []
    y = ya
    y1 = y1a
    for i in range(n):
        if i < 2 and fl:
            fl.write(f'x{i + 1} = {x + h}\n'
                     f'z{i + 1} = z{i} + hf(x{i},y{i}, y\'{i}) = '
                     f'{y1} + {h} * ({np.sin(3 * x)} - {y}) = {y1 + h * f(x, y, y1)}\n'
                     f'y{i + 1} = y{i} + hz = {y + h * (y1 + h * f(x, y, y1))}\n\n')
        deltaYk.append(h * f(x, y, y1))
        y1 += h * f(x, y, y1)
        y += h * y1
        x += h
        x_res.append(x)
        y_res.append(y)
    deltaYk.append('-')
    x_res = [round(x, 5) for x in x_res]
    return x_res, y_res, deltaYk


def runge_kutta(f, xa, xb, ya, y1a, h, fi=None):
    n = int((xb - xa) / h)
    x = xa
    y = ya
    z = y1a
    x_res = [x]
    y_res = [y]
    z_res = [z]
    deltaY = []
    deltaZ = []
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
        deltaY.append((k1 + 2 * k2 + 2 * k3 + k4) / 6)
        deltaZ.append((l1 + 2 * l2 + 2 * l3 + l4) / 6)

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        z += (l1 + 2 * l2 + 2 * l3 + l4) / 6
        if fi and i == 1:
            fi.write(f'\nRunge-Kutt:\n\n')
            fi.write(f"K_1 = {k1}\n"
                     f"L_1 = {l1}\n"
                     f"K_2 = {k2}\n"
                     f"L_2 = {l2}\n"
                     f"K_3 = {k3}\n"
                     f"L_3 = {l3}\n"
                     f"K_4 = {k4}\n"
                     f"L_4 = {l4}\n\n")
            fi.write(f"deltaY_0 = {deltaY[0]}\n"
                     f"deltaZ_0 = {deltaZ[0]}\n\n")
            fi.write(f"y1 = {y}\n")
            fi.write(f"z1 = {z}\n\n")
        x_res.append(x)
        y_res.append(y)
        z_res.append(z)
    deltaY.append(0)
    deltaZ.append(0)
    x_res = [round(x, 5) for x in x_res]
    return x_res, y_res, z_res, deltaY, deltaZ


def adams(f, x, y, h, n, z):
    z = z[:4] + [0] * (len(z) - 4)
    for i in range(3, n):
        z[i + 1] = z[i] + h / 24 * (55 * f(x[i], y[i], z[i]) - 59 * f(x[i - 1], y[i - 1], z[i - 1]) +
                                    37 * f(x[i - 2], y[i - 2], z[i - 2]) - 9 * f(x[i - 3], y[i - 3], z[i - 3]))
        tmp = y[i] + h / 24 * (55 * z[i] - 59 * z[i - 1] + 37 * z[i - 2] - 9 * z[i - 3])
        x.append(x[-1] + h)
        y.append(tmp)
    x = [round(i, 2) for i in x]
    return x, y, z


def print_result_table(f, name, res, runge_y, p):
    if 'Euler' == name:
        table = PrettyTable(['x', 'y', 'delta_Yk', 'y_src', 'eps_k',
                             'Romberg\'s method error'])
        for x, y, delta, yr in zip(*res, runge_y):
            tmp = orig_foo(x)
            table.add_row([round(x, 8), round(y, 8), delta, round(tmp, 8), abs(y - tmp),
                           abs(y + abs(y - yr) / (0.5 ** p - 1) - tmp)])
        f.write(f'\n{str(table)}\n')
    elif 'Runge-Kutta' == name:

        table = PrettyTable(['x', 'y', 'z', 'delta_Yk', 'delta_Zk', 'y_src', 'eps_k',
                             'Romberg\'s method error'])
        for x, y, z, deltaY, deltaZ, yr in zip(*res, runge_y):
            tmp = orig_foo(x)
            table.add_row([round(x, 5), round(y, 5), round(z, 5), round(deltaY, 5),
                           round(deltaZ, 5), round(tmp, 5), abs(round(y - tmp, 6)),
                           abs(y + abs(y - yr) / (0.5 ** p - 1) - tmp)])
        f.write(f'\n{str(table)}\n')
    elif 'Adams' == name:
        f.write(f'\nAdams:\n')
        table = PrettyTable(['x_k', 'y_k', 'f(x_k, y_k)', 'y_src', 'eps_k',
                             'Romberg\'s method error'])
        for x, y, z, yr in zip(*res, runge_y):
            tmp = orig_foo(x)
            table.add_row([round(x, 5), round(y, 5), round(z, 5), round(tmp, 5), abs(round(y - tmp, 6)),
                           abs(round(y + abs(y - yr) / (0.5 ** p - 1) - tmp, 7))])
        f.write(f'\n{str(table)}\n')


def main():

    a = 0
    b = 1
    y0 = 1
    y10 = 1
    step = 0.1
    f = open("output", 'w')

    res1 = euler(foo, a, b, y0, y10, step, f)
    res1_half_h = euler(foo, a, b, y0, y10, step / 2)
    h_half = [y for x, y in zip(res1_half_h[0], res1_half_h[1]) if x in res1[0]]
    print_result_table(f, 'Euler', res1, h_half, 1)

    res2 = runge_kutta(foo, a, b, y0, y10, step, f)
    z = res2[2]
    res2_half_h = runge_kutta(foo, a, b, y0, y10, step / 2)
    z_half_h = res2_half_h[2]
    h_half = [y for x, y in zip(res2_half_h[0], res2_half_h[1]) if x in res2[0]]
    print_result_table(f, 'Runge-Kutta', res2, h_half, 4)

    res3 = adams(foo, res2[0][:4], res2[1][:4], step, int((b - a) / step), z)
    res3_half_h = adams(foo, res2[0][:4], res2[1][:4], step / 2,
                        int((b - a) / (step / 2)), z_half_h)
    h_half = [y for x, y in zip(res3_half_h[0], res3_half_h[1]) if x in res3[0]]
    print_result_table(f, 'Adams', res3, h_half, 4)

    f.close()


if __name__ == "__main__":
    main()
