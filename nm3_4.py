"""
Output:
0.433015
-0.13397000000000003
"""

def find_start(x, p):
    for i in range(0, len(p) - 1):
        if p[i] <= x and x <= p[i + 1]:
            return i

def df1(x, y, x0):
    i = find_start(x0, x)
    elem1 =  (y[i + 1] - y[i]) / (x[i + 1] - x[i])
    elem2 =  ((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) - elem1) / (x[i + 2] - x[i]) * (2 * x0 - x[i] - x[i + 1])
    return elem1 + elem2

def df2(x, y, x0):
    i = find_start(x0, x)
    elem1 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1])
    elem2 = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
    return 2 * (elem1 - elem2) / (x[i + 2] - x[i])


x = [-1.0,0.0,1.0,2.0,3.0]
y = [-0.5,0.0,0.5,0.86603,1.0]

x_star = 1.0

print(df1(x, y, x_star))
print(df2(x, y, x_star))