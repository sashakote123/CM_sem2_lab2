import numpy as np
from tabulate import tabulate
import math


def setMatrix2(n, m):
    h = 1 / n
    k = 1 / m
    a = 1 / h ** 2
    b = 1 / k ** 2
    A = -2 * (a + b)

    matrix = []
    for i in range((m - 1) * (n - 1)):
        matrix.append([0] * (m - 1) * (n - 1))

    for i in range((n - 1) * (m - 1)):
        matrix[i][i] = A
    for i in range((n - 1) * (m - 2)):
        matrix[i][i + n - 1] = b
    for i in range((n - 1) * (m - 2)):
        matrix[i + n - 1][i] = b
    for i in range((n - 1) * (m - 1) - 1):
        matrix[i + 1][i] = a
    for i in range((n - 1) * (m - 1) - 1):
        matrix[i][i + 1] = a

    for i in range(1, (m - 1)):
        matrix[i * (n - 1) - 1][i * (n - 1)] = 0

    for i in range(1, (m - 1)):
        matrix[i * (n - 1)][i * (n - 1) - 1] = 0

    return matrix


def set_solve_vector(n, m):
    h = 1 / n
    k = 1 / m

    m1_list = []
    m2_list = []
    m3_list = []
    m4_list = []
    f_list = []

    y = 0 + k
    x = -1 + h
    for i in range(m - 1):
        m1_list.append(np.exp(-y))
        m2_list.append(1)
        y += k

    for i in range(n - 1):
        m3_list.append(1)
        m4_list.append(np.exp(x))
        x += h

    x = -1
    y = 0
    for i in range(m - 1):
        x += h
        y = 0
        for j in range(n - 1):
            y += k
            f_list.append(np.cosh(x ** 2 * y))

    vector = []
    for i in range((n - 1) * (m - 1)):
        vector.append(-f_list[i])

    for i in range(len(m3_list)):
        vector[i] -= m3_list[i] / k ** 2
        vector[-1 - i] -= m4_list[-1 - i] / k ** 2

    for i in range(len(m1_list)):
        vector[i * len(m3_list)] -= m1_list[i] / h ** 2
        vector[i * len(m3_list) + len(m3_list) - 1] -= m2_list[i] / h ** 2

    return vector

def set_solve_vector3(n, m):
    h = 1 / n
    k = 1 / m

    m1_list = []
    m2_list = []
    m3_list = []
    m4_list = []
    f_list = []

    y = 0 + k
    x = -1 + h
    for i in range(m - 1):
        m1_list.append(np.sin(np.pi * y))
        m2_list.append(abs(np.sin(2 * np.pi * y)))
        y += k

    for i in range(n - 1):
        m3_list.append(-x * (x + 1))
        m4_list.append(-x * (x + 1))
        x += h

    x = -1
    y = 0
    for i in range(m - 1):
        x += h
        y = 0
        for j in range(n - 1):
            y += k
            f_list.append(np.cosh(x ** 2 * y))

    vector = []
    for i in range((n - 1) * (m - 1)):
        vector.append(-f_list[i])

    for i in range(len(m3_list)):
        vector[i] -= m3_list[i] / k ** 2
        vector[-1 - i] -= m4_list[-1 - i] / k ** 2

    for i in range(len(m1_list)):
        vector[i * len(m3_list)] -= m1_list[i] / h ** 2
        vector[i * len(m3_list) + len(m3_list) - 1] -= m2_list[i] / h ** 2

    return vector


def linear_interpolation(solve_vector, h):
    x = -1 + h
    interp = []
    for j in range(int(np.sqrt(len(solve_vector)))):
        for i in range(int(np.sqrt(len(solve_vector)))):
            interp.append(-x * (x + 1))
            x += h
        x = -1 + h
    return interp


def Zeydel_solve_test(n, m, N_max, eps, omega):
    # N_max = 10000
    S = 0
    # eps = 1e-12

    a = -1
    b = 0
    c = 0
    d = 1
    v_new = 0
    exit = False
    h2 = -(n / (b - a)) ** 2
    k2 = -(m / (d - c)) ** 2
    a2 = -2 * (h2 + k2)

    f_list = []
    f_test_list = []
    v_list = []
    r_list = []

    h = 1 / n
    k = 1 / m

    for i in range((m + 1)):
        f_list.append([0] * (n + 1))
        v_list.append([0] * (n + 1))
        r_list.append([0] * (n + 1))
        f_test_list.append([0] * (n + 1))

    # v_list.reverse()
    # f_list.reverse()

    y = 0
    for i in range(m + 1):
        x = -1
        for j in range(n + 1):
            f_list[i][j] = (np.exp(x * y) * x ** 2 + np.exp(x * y) * y ** 2)
            v_list[i][j] = -np.exp(x * y)
            r_list[i][j] = -np.exp(x * y)
            f_test_list[i][j] = -np.exp(x * y)
            x += h
        y += k

    for i in range(1, m):
        for j in range(1, n):
            v_list[i][j] = 0
    # f_list.reverse()
    # print("f list:", f_list)
    # print("v list:", v_list)
    while exit != True:
        eps_max = 0
        for j in range(1, n):
            for i in range(1, m):
                v_old = v_list[i][j]
                v_new = -omega * (
                        h2 * (v_list[i + 1][j] + v_list[i - 1][j]) + k2 * (v_list[i][j + 1] + v_list[i][j - 1]))
                v_new = v_new + (1 - omega) * a2 * v_list[i][j] + omega * f_list[i][j]
                v_new /= a2
                eps_cur = abs(v_old - v_new)
                if eps_cur > eps_max:
                    eps_max = eps_cur
                v_list[i][j] = v_new
        S += 1
        if (eps_max < eps) or (S >= N_max):
            exit = True

    for i in range(m + 1):
        for j in range(n + 1):
            r_list[i][j] = abs(r_list[i][j] - v_list[i][j])
    matrix = np.array(r_list)
    max_index = np.argmax(r_list, axis=None)
    row_index, col_index = np.unravel_index(max_index, matrix.shape)

    Xmax = -1 + col_index * h
    Ymax = 0 + row_index * k

    return v_list, max(max(r_list)), S, eps_max, f_list, f_test_list, Xmax, Ymax, r_list


def Zeydel_solve(n, m, N_max, eps, omega):
    # N_max = 10000
    S = 0
    # eps = 1e-12

    a = -1
    b = 0
    c = 0
    d = 1
    v_new = 0
    exit = False
    h2 = -(n / (b - a)) ** 2
    k2 = -(m / (d - c)) ** 2
    a2 = -2 * (h2 + k2)

    f_list = []
    v_list = []
    r_list = []

    h = 1 / n
    k = 1 / m

    for i in range((m + 1)):
        f_list.append([0] * (n + 1))
        v_list.append([0] * (n + 1))
        r_list.append([0] * (n + 1))

    # v_list.reverse()
    # f_list.reverse()

    # y = c + k
    # for i in range(1, m):
    #    x = -1 + h
    #    for j in range(1, n):
    #        f_list[i][j] = -np.cosh(x ** 2 * y)
    #        x += h
    #    y += k

    y = 0
    for i in range(m + 1):
        x = -1
        for j in range(n + 1):
            f_list[i][j] = np.cosh(x ** 2 * y)
            x += h
        y += k

    y = c
    for i in range(n + 1):
        v_list[0][i] = np.sin(np.pi * y)
        v_list[-1][i] = abs(np.sin(2 * np.pi * y))
        y += k

    x = -1
    for i in range(m + 1):
        v_list[i][0] = -x * (x + 1)
        v_list[i][-1] = -x * (x + 1)
        x += h

    # f_list.reverse()
    # print("f list:", f_list)
    # print("v list:", v_list)
    while exit != True:
        eps_max = 0
        for j in range(1, n):
            for i in range(1, m):
                v_old = v_list[i][j]
                v_new = -omega * (
                        h2 * (v_list[i + 1][j] + v_list[i - 1][j]) + k2 * (v_list[i][j + 1] + v_list[i][j - 1]))
                v_new = v_new + (1 - omega) * a2 * v_list[i][j] + omega * f_list[i][j]
                v_new /= a2
                eps_cur = abs(v_old - v_new)
                if eps_cur > eps_max:
                    eps_max = eps_cur
                v_list[i][j] = v_new
        S += 1
        if (eps_max < eps) or (S >= N_max):
            exit = True

    return v_list, r_list, S, eps_max, f_list


def splitMatrix(n, m):
    matrix = setMatrix2(n, m)
    MatrixL = []
    MatrixR = []
    MatrixD = []
    MatrixD1 = []
    for i in range((m - 1) * (n - 1)):
        MatrixL.append([0] * (m - 1) * (n - 1))
        MatrixR.append([0] * (m - 1) * (n - 1))
        MatrixD.append([0] * (m - 1) * (n - 1))
        MatrixD1.append([0] * (m - 1) * (n - 1))

    for i in range((m - 1) * (n - 1)):
        for j in range((m - 1) * (n - 1)):
            if i < j:
                MatrixR[i][j] = matrix[i][j]
            if i > j:
                MatrixL[i][j] = matrix[i][j]
        MatrixD[i][i] = matrix[i][i]
        MatrixD1[i][i] = 1 / matrix[i][i]
    return MatrixD, MatrixL, MatrixR, MatrixD1


def getB(n, m):
    D, L, R, D1 = splitMatrix(n, m)

    tmp = []
    for i in range((m - 1) * (n - 1)):
        tmp.append([0] * (m - 1) * (n - 1))
    for i in range((m - 1) * (n - 1)):
        for j in range((m - 1) * (n - 1)):
            tmp[i][j] = L[i][j] - R[i][j]
    x = np.dot(D1, tmp)
    return np.linalg.eig(x)[0]


def getOmega(n, m):
    v = getB(n, m)
    for i in range((n - 1) * (m - 1)):
        v[i] = abs(v[i])
    return 2 / (1 + np.sqrt(1 - (max(v) ** 2).real))

n = 5
m = 6
matrix = setMatrix2(n, m)
solve = set_solve_vector3(n, m)
v = Zeydel_solve(n, m, 1000, 1e-12, 1.5)[0]
v3 = np.linalg.solve(matrix,solve)

def nevyazka(matrix, solve, v):
    v2 = []
    v3 = np.linalg.solve(matrix, solve)
    for i in range(1, len(v[0]) - 1):
        for j in range(1, len(v) - 1):
            v2.append(v[j][i])
    dot = np.dot(matrix, v2)
    nev = []
    for i in range(len(solve)):
        nev.append(abs(solve[i] - dot[i]))
    return nev


#print(max(nevyazka(matrix, solve, v)))

