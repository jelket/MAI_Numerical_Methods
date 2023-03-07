"""
Output:
Method of rectangles = -0.05450105450105448	Step = 0.5
Method of rectangles = -0.05794799368631646	Step = 0.25
Error: 0.00002532658526288589
Romberg's method: -0.05909697341473712

Method of trapeziums = -0.06845238095238095	Step = 0.5
Method of trapeziums = -0.061476717726717735	Step = 0.25
Error: 0.00002919665149666434
Romberg's method: -0.05915149665149667

Simpson's method  = -0.05952380952380953	Step = 0.5
Simpson's method  = -0.05915149665149664	Step = 0.25
Error: 0.00000437579334244126
Romberg's method: -0.059126675793342444
"""

import numpy as np


def f(x):
    return x / (2 * x + 5)


def simpson(a, b, h, n):
    s = f(a) + f(b)
    for i in range(1, n, 2):
        s += 4 * f(a + i * h)
    for i in range(2, n - 1, 2):
        s += 2 * f(a + i * h)
    return s * h / 3


def main():

    a, b, h1, h2 = -1, 1, 0.5, 0.25

    x1 = np.linspace(a, b, int((b - a) / h1 + 1))
    x2 = np.linspace(a, b, int((b - a) / h2 + 1))

    y_trap1 = [f(i) for i in x1]
    y_trap2 = [f(i) for i in x2]

    rect1 = h1 * sum([f((i + j) / 2) for i, j in zip(x1, x1[1:])])
    rect2 = h2 * sum([f((i + j) / 2) for i, j in zip(x2, x2[1:])])

    trap1 = h1 * (y_trap1[0] / 2 + sum(y_trap1[1:-1]) + y_trap1[-1] / 2)
    trap2 = h2 * (y_trap2[0] / 2 + sum(y_trap2[1:-1]) + y_trap2[-1] / 2)

    simps1 = simpson(a, b, h1, int((b - a) / h1))
    simps2 = simpson(a, b, h2, int((b - a) / h2))

    true_val = -0.0591223
    print(f'Method of rectangles = {rect1}\tStep = {h1}')
    print(f'Method of rectangles = {rect2}\tStep = {h2}')
    print(f'Error: {abs(rect1 + ((rect1 - rect2) / (0.5**2 - 1)) - true_val):.20f}')
    print(f'Romberg\'s method: {rect1 + ((rect1 - rect2) / (0.5**2 - 1))}')
    print()
    print(f'Method of trapeziums = {trap1}\tStep = {h1}')
    print(f'Method of trapeziums = {trap2}\tStep = {h2}')
    print(f'Error: {abs(trap1 + ((trap1 - trap2) / (0.5**2 - 1)) - true_val):.20f}')
    print(f'Romberg\'s method: {trap1 + ((trap1 - trap2) / (0.5**2 - 1))}')
    print()
    print(f"Simpson's method  = {simps1}\tStep = {h1}")
    print(f"Simpson's method  = {simps2}\tStep = {h2}")
    print(f'Error: {abs(simps1 + ((simps1 - simps2) / (0.5**4 - 1)) - true_val):.20f}')
    print(f'Romberg\'s method: {simps1 + ((simps1 - simps2) / (0.5**4 - 1))}')


if __name__ == "__main__":
    main()
