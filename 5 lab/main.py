import numpy as np
import scipy.linalg as sp

# constants
k = 23
m = 6
eps = 0.0001
a = np.array([
            [12+k, 2, m/4, 1, 2],
            [4, 113+k, 1, m/10, m-4],
            [1, 2, -24-k, 3, 4],
            [1, 2/m, 4, 33+k, 4],
            [-1, 2, -3, 3+m, -44-k]
            ], float)
b = np.array([1, 2, 3, 4, 5], float)


def is_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol, atol)


def sqrroot(A, B):

    if is_symmetric(A):
        t, tlow = sp.cho_factor(A)
        X = sp.cho_solve((t, tlow), B)
    else:
        print("матрица не семетрична")
        return None

    return X


def yakobi(A, b, tolerance=eps, max_iterations=1000):

    x = np.zeros_like(b, np.double)
    D = np.diag(A)
    T = A - np.diagflat(D)
    countiter = 0

    for k in range(max_iterations):
        countiter += 1
        x_old = x.copy()
        x[:] = (b - np.dot(T, x)) / D
        if np.linalg.norm(x - x_old) / np.linalg.norm(x) < tolerance:
            break

    return x, countiter


if __name__ == '__main__':
    print("Матрица A")
    print(a)
    print("Метод квадратного корня:", end=' ')
    resSqrroot = sqrroot(a, b)
    if not resSqrroot is None:
        print("X:")
        print(resSqrroot)
        print("A*X=B верно?")
        print(np.allclose(np.dot(a, resSqrroot), b))
    else:
        pass
    print("Метод Якоби:")
    resYakobi, iters = yakobi(a, b)
    print("X: ", resYakobi)
    print("Колличестов шагов: ", iters)
    #print(iters)
    print("A*X=B Верно?", (np.allclose(np.dot(a, resYakobi), b)))
    slv = np.linalg.solve(a, b)
    print("Метод квадратного корня:", slv)
    print("A*X=B с помощью метода квадратного корня верно?", (np.allclose(np.dot(a, slv), b)))