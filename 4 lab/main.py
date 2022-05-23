import numpy as np
import scipy.optimize as sp

N = len('IvanovRuslanAlbertovich')  # 23
m = len('Ruslan')  # 6
variant = N % 4  # 23%4 = 3
a1 = 3 + 0.1 * m  # a1=3.6
b1 = 0.4 + 0.03 * m  # b1=0.58
a2 = 1 + 0.5 * m  # a2=4
eps1 = 10 ** (-4)  #
eps2 = 10 ** (-2)  #


def func1(x, arg1, arg2):
    return (x ** 5) - (arg1 * x) + arg2  # x**5 - a*x +b


def sys1(params):
    x, y = params
    eq1 = a2 * x ** 3 - y ** 2 - 1  # ax^3 - y^2 - 1 = 0
    eq2 = x * y ** 3 - y - 4  # xy^3 - y - 4 =0
    return eq1, eq2  #


if __name__ == '__main__':
    print("Первое задание: ")
    result = sp.fsolve(func1, x0=np.array([1]), args=(a1, b1), maxfev=100000, xtol=eps1)
    print('Найдена фиксированная точка: ', end="")
    print(result)
    print('Значение функции с фиксированной точкой: ', end="")
    print(func1(result, a1, b1))

    print("Второе задание: ")
    result2 = sp.fsolve(sys1, x0=np.array([0.0, 1.0]), xtol=eps2, maxfev=1000)
    xres, yres = result2
    print("x = ", end=" ")
    print(xres)
    print("y = ", end=" ")
    print(yres)
    print(sys1((xres, yres)))
