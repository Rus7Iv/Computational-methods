import numpy as np


def horner(poly, n, x):
    result = poly[0]
    for i in range(1, n):
        result = result * x + poly[i]

    return result


def main():
    poly = [1.234, 0.387, 1.4789, 1.0098, 1.222, -2.345]
    print("\nИсходный многочлен: P(x) = 1.234x^5 + 0.387x^4 + 1.4789x^3 + 1.0098x^2 + 1.222x-2.345")
    x = 0.23426
    n = len(poly)

    print("Найденные значения многочлена P(x) в точке x = 0.23426:")
    print(4 * " ", "- по схеме Горнера: ", horner(poly, n, x))
    print(4 * " ", "- с помощью библиотеки Numpy: ", np.polyval(poly, x))
    c = 0.987
    delta = x - c

    print(" Коэффициенты многочлена, после деления на x - c = ", delta)
    for i in range(0, n):
        print(5 * " ", "x ^", n - i - 1, ":", poly[i] / delta)


main()
