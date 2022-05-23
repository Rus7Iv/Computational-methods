import numpy as np
from numpy import matrix

k = 28
m = 6
A = np.array(
    [[12 + k, 2, m / 4, 1, 2],
     [4, 113 + k, 1, m / 10, m - 4],
     [1, 2, -24 - k, 3, 4],
     [1, 2 / m, 4, 33 + k, 4],
     [-1, 2, -3, 3 + m, -44 - k]])
B = np.array([1, 2, 3, 4, 5])


def decomposeToLinearEquations(a):
    lu_matrix = np.matrix(np.zeros([a.shape[0], a.shape[1]]))
    n = a.shape[0]

    for k in range(n):
        for j in range(k, n):
            lu_matrix[k, j] = a[k, j] - lu_matrix[k, :k] * lu_matrix[:k, j]
        for i in range(k + 1, n):
            lu_matrix[i, k] = (a[i, k] - lu_matrix[i, :k] *
                               lu_matrix[:k, k]) / lu_matrix[k, k]

    return lu_matrix


def getLowerTriangleMatrix(m):
    L = m.copy()
    for i in range(L.shape[0]):
        L[i, i] = 1
        L[i, i + 1:] = 0
    return np.matrix(L)


def getUpperTriangleMatrix(m):
    U = m.copy()
    for i in range(1, U.shape[0]):
        U[i, :i] = 0
    return U


def calculateY(lu_matrix, b):
    y = np.matrix(np.zeros([lu_matrix.shape[0], 1]))
    for i in range(y.shape[0]):
        y[i] = b[i] - lu_matrix[i, :i] * y[:i]
    return y


def calculateX(lu_matrix, y):
    x = np.matrix(np.zeros([lu_matrix.shape[0], 1]))
    for i in range(1, x.shape[0] + 1):
        x[-i] = (y[-i] - lu_matrix[-i, -i:] * x[-i:]) / lu_matrix[-i, -i]
    return x


def kholetsky(A, B):
    decomposedMatrix = decomposeToLinearEquations(A)

    print(" МЕТОД ХОЛЕЦКОГО: \n")
    print(" Исходная матрица: \n")
    print(matrix(A), "\n")
    print(" Решение: \n")
    print(" Полученная верхнетреугольная матрица: \n")
    print(matrix(getUpperTriangleMatrix(decomposeToLinearEquations(A))), "\n")
    print(" Полученная нижнетреугольная матрица: \n")
    print(matrix(getLowerTriangleMatrix(decomposeToLinearEquations(A))), "\n")

    print(" Вычисленные значения y: \n")
    print(calculateY(decomposedMatrix, B), "\n")

    print(" Вычисленные значения x: \n")
    print(calculateX(decomposedMatrix, calculateY(decomposedMatrix, B)), "\n")


def zeidel(A, B, eps):
    n = len(A)
    x = np.zeros(n)

    converge = False
    while not converge:
        x_new = np.copy(x)
        for i in range(n):
            s1 = sum(A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (B[i] - s1 - s2) / A[i][i]

        converge = np.sqrt(
            sum((x_new[i] - x[i]) ** 2 for i in range(n))) <= eps
        x = x_new

    print(" МЕТОД ЗЕЙДЕЛЯ: \n")
    print(" Исходная матрица: \n")
    print(matrix(A), "\n")
    print(" Вычисленные значения x: \n")
    print(matrix(x), "\n")


def main():
    kholetsky(A, B)
    zeidel(A, B, eps=0.01)
    print(" Вычисленные с помощью библиотеки numpy значения x:")
    print(matrix(np.linalg.solve(A, B)), "\n")


main()
