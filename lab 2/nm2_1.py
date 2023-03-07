"""
Output:
q = 0.6301279573390394
lambda = -0.8132227352659079

phi' < 1
q < 1

Iter:
x0 = 1.25
x1 = 1.5069086342344487, eps = 0.4376792356304313
x2 = 1.5648348676802541, eps = 0.09868531531864919
x3 = 1.5727460695657183, eps = 0.013477821811890169
x4 = 1.5737005477535821, eps = 0.0016260850279907106
x5 = 1.573813732668428, eps = 0.000192826088396256
Newton's:
x5 = 1.25
x1 = 1.6210487834750746, eps = 0.37104878347507464
x2 = 1.574369411050435, eps = 0.04667937242463971
x3 = 1.5738290003277182, eps = 0.000540410722716711


Solution:
EPS = 0.001
Iter: x = {self.iter_x}
Steps = 5
Newton: x = 1.5738290003277182
Steps = 3

Process finished with exit code 0

"""

import numpy as np


class Solver:
    def __init__(self, eps):
        self.eps = eps
        self.area = (1, 3)
        self.x0 = 1.25
        self.lmbd = self.calc_lambda()
        self.q = self.calc_q()
        self.k_iter = 0
        self.k_newton = 0
        self.check()
        self.iter_x = self.iter_method()
        self.newtons_x = self.newtons_method()

    def check(self):
        print(f'q = {self.q}')
        print(f'lambda = {self.lmbd}\n')
        x = np.linspace(self.area[0], self.area[1], 10000)
        y = [self.phi_derivative(i) for i in x]
        if all([i < 1 for i in y]):
            print('phi\' < 1')
        else:
            print('phi\' >= 1')
        if self.q < 1:
            print('q < 1\n')
        else:
            print('q >= 1\n')

    @staticmethod
    def f(x):
        return 2**x - x**2 - 0.5

    def phi(self, x):
        return x - self.lmbd * self.f(x)

    def phi_derivative(self, x):
        return 1 - self.lmbd * self.f_derivative(x)

    @staticmethod
    def f_derivative(x):
        return 2**x * np.log(2) - 2*x

    @staticmethod
    def f_2derivative(x):
        return 2**x * np.log(2)**2 + 2**x * 0.5 - 2

    def calc_q(self):
        x = np.linspace(self.area[0], self.area[1], 10000)
        y = [abs(self.phi_derivative(i)) for i in x]
        q = np.max(y)
        return q

    def calc_lambda(self):
        flag = None
        x = np.linspace(self.area[0], self.area[1], 10000)
        y = [self.f_derivative(i) for i in x]

        if all([np.sign(i) == -1 for i in y]):
            flag = -1
        elif all([np.sign(i) == 1 for i in y]):
            flag = 1
        else:
            print('Error: Derivative change sign\n')
            exit(-1)

        y = [abs(self.f_derivative(i)) for i in x]
        return flag / np.max(y)

    def iter_method(self):
        x_old = self.x0
        print(f'Iter:\nx{self.k_iter} = {x_old}')
        while True:
            self.k_iter += 1
            x_new = self.phi(x_old)
            print(f'x{self.k_iter} = {x_new}, eps = {abs(x_new - x_old) * self.q / (1 - self.q)}')
            if abs(x_new - x_old) * self.q / (1 - self.q) < self.eps:
                return x_new
            else:
                x_old = x_new

    def newtons_method(self):
        x_old = self.x0
        print(f'Newton\'s:\nx{self.k_iter} = {x_old}')
        while True:
            self.k_newton += 1
            x_new = x_old - self.f(x_old) / self.f_derivative(x_old)
            print(f'x{self.k_newton} = {x_new}, eps = {abs(x_new - x_old)}')
            if abs(x_new - x_old) < self.eps:
                return x_new
            else:
                x_old = x_new

    def print_solution(self):
        print("\n")
        print("Solution:")
        print(f'EPS = {self.eps}')
        print('Iter: x = {self.iter_x}')
        print(f'Steps = {self.k_iter}')
        print(f'Newton: x = {self.newtons_x}')
        print(f'Steps = {self.k_newton}')


def main():

    sol = Solver(0.001)
    sol.print_solution()


if __name__ == "__main__":
    main()
