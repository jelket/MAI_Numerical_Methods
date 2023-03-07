"""
Output:
Simple iterations:
[7.999986531267736, 4.000008718066204, 3.0000166419957215, 8.99997759770737] 36
Zeidel:
[7.999958874657556, 3.99997127811487, 2.9999759805373163, 8.999974804931705] 16
"""


import copy
import math


def simple(A, b, err):
    x0 = [None] * len(A)
    x = [0] * len(A)
    num_of_it = 0
    while True:
        for i in range(len(A)):
            s = 0
            for j in range(len(A)):
                if i != j:
                    s += A[i][j] * x[j]
            x0[i] = (b[i] - s) / A[i][i]
        num_of_it += 1
        if error(x,x0,err):
            break
        x = copy.copy(x0)
    print("Simple iterations:")
    print(x0, num_of_it)


def zeidel(A, b, err):
    x0 = [None] * len(A)
    x = [0] * len(A)
    num_of_it = 0
    while True:
        for i in range(len(A)):
            s = 0
            for j in range(len(A)):
                if j < i:
                    s += A[i][j] * x0[j]
                elif i != j:
                    s += A[i][j] * x[j]
            x0[i] = (b[i] - s) / A[i][i]
        num_of_it += 1
        if error(x,x0,err):
            break
        x = copy.copy(x0)
    print("Zeidel:")
    print(x0, num_of_it)


def error(x,x0, err):
    if x0[0] == None:
        return False
    res = [0] * len(x)
    for i in range(len(x)):
        res[i] = pow(x[i] - x0[i],2)
    if math.sqrt(sum(res)) > err:
        return False
    else:
        return True


err = 10**-4
A = [[19, -4, -9, -1],
         [-2, 20, -2, -7],
         [6, -5, -25, 9],
         [0, -3, -9, 12]]
b = [100, -5, 34, 69]

simple(A, b, err)
zeidel(A, b, err)
