import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from backend import *

n = 16
m = 16

v, r, S, eps, f_list, f, Xmax, Ymax, r_list = Zeydel_solve_test(n, m, 1000, 1e-12, 1.5)



v2 = Zeydel_solve(n, m, 100, 1e-12, 1.5)[0]
v3 = Zeydel_solve(2 * n, 2 * m, 100, 1e-12, 1.5)[0]

r_list = []
for i in range((m + 1)):
    r_list.append([0] * (n + 1))

for i in range(m + 1):
    for j in range(n + 1):
        r_list[i][j] = abs(v2[i][j] - v3[(2 * i)][(2 * j)])

print(v)


def getCoordinates(n, m, v):
    h = 1 / n
    k = 1 / m

    X = []
    Y = []
    Z = []
    for i in range(len(v)):
        for j in range(len(v[i])):
            X.append(-1 + i * h)
            Y.append(j * k)
            Z.append(v[i][j])
    return X, Y, Z


def plotGraph(X,Y,Z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(X, Y, Z)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()


