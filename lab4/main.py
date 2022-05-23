from scipy.optimize import fsolve

def newton(f, Df, x0, epsilon, max_iter):
    xn = x0

    for n in range(0, max_iter):
        fxn = f(xn)
        if abs(fxn) < epsilon:
            return xn
        Dfxn = Df(xn)
        if Dfxn == 0:
            return None
        xn = xn - fxn / Dfxn
    
    return None

def equations(vars):
    x, y = vars
    eq1 = 1.3 * x ** 3 - y ** 2 - 1
    eq2 = x * y ** 3 - y - 4
    
    return [eq1, eq2]

def main():
    m = 6
    a = 3 + 0.1 * m
    b = 0.4 + 0.03 * m
    eps = 0.0001

    f = lambda x: x ** 5 - a * x + b
    Df = lambda x: 5 * x ** 4 - a

    print('\nЗАДАНИЕ А:')
    print(' ИСХОДНОЕ УРАВНЕНИЕ: x^5 - ax + b = 0, ГДЕ ')
    print('     a = ', a)
    print('     b = ', b)
    print('\n ПОЛОЖИТЕЛЬНЫЕ КОРНИ УРАВНЕНИЯ, НАЙДЕННЫЕ С ПОМОЩЬЮ МЕТОДА НЬЮТОНА:')

    result = newton(f, Df, 0, eps, 10)
    print('     x1 = ', result)
    
    result = newton(f, Df, 1, eps, 10)
    print('     x2 = ', result)

    print('\nЗАДАНИЕ Б:')
    print(' ИСХОДНАЯ СИСТЕМА УРАВНЕНИЙ:')
    print('     ax^3 - y^2 - 1 = 0 ')
    print('     xy^3 - y - 4 = 0, ')
    print(' ГДЕ a = 1.3')
    print('\n РЕШЕНИЕ СИСТЕМЫ УРАВНЕНИЙ: ')

    x, y =  fsolve(equations, (0, 2))
    
    print('     x = ', x)
    print('     y = ', y, '\n')


main()
