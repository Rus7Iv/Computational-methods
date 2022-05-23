from scipy import integrate
import numpy as np

a = 0
b = np.pi
eps = 10 ** (-4)


def f(x):
    return np.sqrt(1 + np.sin(1 / (1 + x ** 2)) ** 2)


def simpson(f, a, b, eps):
    result = ((b - a) / 6) * (f(a) + 4 * f((a+b) / 2) + f(b))
    x = a + eps

    while x <= b - eps:
        result += f(x)
        x += eps

    return eps * result


def main():
    print("\nДано: [a,b] = [0,pi]")
    print("     f(x): sqrt(1 + sin^2(1/(1+x^2)) \n")
    print("РЕЗУЛЬТАТ ВЫЧИСЛЕНИЯ ИНТЕГРАЛОВ: ")
    print(" 1. По формуле Симпсона: " +
          str(simpson(f, a, b, eps)))
    print(" 2. По квадратурной формуле Гаусса с 3-мя узлами: " +
          str(integrate.fixed_quad(f, a, b, n=3)[0]) + "\n")


main()
