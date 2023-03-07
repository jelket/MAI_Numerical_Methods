"""
Output:
Iter:
x0 = [2.5, 0.5]
x1 = [3.04887218 0.64661654], eps = 2.461437097899386
x2 = [2.93141567 0.6170041 ], eps = 0.5267379512665964
x3 = [2.97623247 0.62850686], eps = 0.20098256531666814
x4 = [2.96037651 0.62444091], eps = 0.0711066457333385
x5 = [2.96616616 0.62593055], eps = 0.02596390155361563
x6 = [2.96407427 0.62539219], eps = 0.009381171453750633
x7 = [2.96483316 0.62558762], eps = 0.0034032981297863448
x8 = [2.96455823 0.62551681], eps = 0.0012329446679830273
x9 = [2.96465789 0.62554248], eps = 0.00044691033381895074
Newton's:
x9 = [2.5, 0.5]
x1 = [3.04887218 0.64661654], eps = 0.5488721804511276
x2 = [2.96657052 0.62610893], eps = 0.08230166419036067
x3 = [2.96463249 0.62553603], eps = 0.0019380266142916547
x4 = [2.96463138 0.62553565], eps = 1.111192574398956e-06

Solution:
q = 0.8176691729323308
Lambda:
[[ 0.03007519  0.30827068]
 [ 0.09022556 -0.07518797]]
EPS = 0.001
Iter: x = [2.96465789 0.62554248]
Steps = 9
Newton: x = [2.96463138 0.62553565]
Steps = 4
"""

import numpy as np
from numpy.linalg import norm, solve, det
from itertools import product


class Solver:
    def __init__(self, eps):
        self.eps = eps
        self.area = ((2.5, 3.25), (0.5, 0.75))
        self.x0 = [2.5, 0.5]
        self.k_iter = 0
        self.k_newton = 0
        self.lmbd = self.calc_lambda()
        self.q = self.calc_q()
        self.iter_x = self.iter_method()
        self.newtons_x = self.newtons_method()

    @staticmethod
    def f1(x1, x2):
        return (x1**2 + 4) * x2 - 8

    @staticmethod
    def f2(x1, x2):
        return (x1 - 1)**2 + (x2 - 1)**2 - 4

    @staticmethod
    def f11(x1, x2):
        return 2 * x1 * x2

    @staticmethod
    def f12(x1):
        return x1**2 + 4

    @staticmethod
    def f21(x1):
        return 2 * x1 - 2

    @staticmethod
    def f22(x2):
        return 2 * x2 - 2

    def phi1(self, x1, x2):
        return x1 - (self.f1(x1, x2) * self.lmbd[0, 0] + self.f2(x1, x2) *
                     self.lmbd[0, 1])

    def phi2(self, x1, x2):
        return x2 - (self.f1(x1, x2) * self.lmbd[1, 0] + self.f2(x1, x2) *
                     self.lmbd[1, 1])

    def phi11(self, x1, x2):
        return 1 - (self.f11(x1, x2) * self.lmbd[0, 0] + self.f21(x1) *
                    self.lmbd[0, 1])

    def phi12(self, x1, x2):
        return -(self.f12(x1) * self.lmbd[0, 0] + self.f22(x2) *
                 self.lmbd[0, 1])

    def phi21(self, x1, x2):
        return -(self.f11(x1, x2) * self.lmbd[1, 0] + self.f21(x1) *
                 self.lmbd[1, 1])

    def phi22(self, x1, x2):
        return 1 - (self.f12(x1) * self.lmbd[1, 0] + self.f22(x2) *
                    self.lmbd[1, 1])

    def phi_derivative(self, x):
        return np.array([[self.phi11(*x), self.phi12(*x)],
                        [self.phi21(*x), self.phi22(*x)]])

    def j(self, x1, x2):
        return [[self.f11(x1, x2), self.f12(x1)], [self.f21(x1), self.f22(x2)]]

    def calc_lambda(self):
        shape = len(self.area)
        current_j = self.j(*self.x0)
        inv_j = np.array([solve(current_j, i) for i in np.eye(shape)])
        return np.transpose(inv_j)

    def calc_q(self):
        x1 = np.linspace(self.area[0][0], self.area[0][1], 100)
        x2 = np.linspace(self.area[1][0], self.area[1][1], 100)
        points = list(product(x1, x2))
        vals = [norm(self.phi_derivative(point), np.inf) for point in points]
        q = np.max(vals)
        return q

    def iter_method(self):
        x_old = self.x0
        print(f'Iter:\nx{self.k_iter} = {x_old}')
        while True:
            self.k_iter += 1
            x_new = np.array([self.phi1(*x_old), self.phi2(*x_old)])
            print(f'x{self.k_iter} = {x_new}, eps = {norm(x_new - x_old, np.inf) * self.q / (1 - self.q)}')
            if norm(x_new - x_old, np.inf) * self.q / (1 - self.q) <= self.eps:
                return x_new
            else:
                x_old = x_new

    def newtons_method(self):
        shape = len(self.area)
        x_old = self.x0
        print(f'Newton\'s:\nx{self.k_iter} = {x_old}')
        while True:
            current_j = self.j(*x_old)
            if det(current_j) == 0:
                print(f'Error: detJ({self.k_newton} == 0)\n')
                exit(-1)
            self.k_newton += 1
            inv_j = np.array([solve(current_j, i) for i in np.eye(shape)])
            x_new = x_old - np.transpose(inv_j) @ np.array([self.f1(*x_old),
                                                            self.f2(*x_old)])
            print(f'x{self.k_newton} = {x_new}, eps = {norm(x_new - x_old, np.inf)}')
            if norm(x_new - x_old, np.inf) <= self.eps:
                return x_new
            else:
                x_old = x_new

    def print_solution(self):
        print()
        print("Solution:")
        print(f'q = {self.q}')
        print(f'Lambda:\n{self.lmbd}')
        print(f'EPS = {self.eps}')
        print(f'Iter: x = {self.iter_x}')
        print(f'Steps = {self.k_iter}')
        print(f'Newton: x = {self.newtons_x}')
        print(f'Steps = {self.k_newton}')


def main():
    sol = Solver(0.001)
    sol.print_solution()


if __name__ == "__main__":
    main()
