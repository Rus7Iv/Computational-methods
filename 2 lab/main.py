import numpy as np
import matplotlib.pyplot as plt

k = 6
m = 6
N = 18

x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
y = [0.2 * N, 0.3 * m, 0.5 * k, 0.6 * N, 0.7 * m, k, 0.8 * N, 1.2 * k, 1.3 * m, N]


def fi1(x):  # Функция, возвращающая значение Фи1
    return 1 - x


def fi2(x, i):
    return x * pow((1 - x), (i - 1))  # Функция, возвращающая значение Фи2


def lagr(x, y, t):  # Формула Лагранжа
    z = 0
    for j in range(len(y)):
        p1 = 1
        p2 = 1
        for i in range(len(x)):
            if i == j:
                p1 = p1 * 1
                p2 = p2 * 1
            else:
                p1 = p1 * (t - x[i])
                p2 = p2 * (x[j] - x[i])
        z = z + y[j] * p1 / p2
    return z


if __name__ == '__main__':
    xnew = np.linspace(np.min(x), np.max(x), 100)  # x
    ynew = [lagr(x, y, i) for i in xnew]  # y

# lagr
    plt.plot(x, y, 'o', xnew, ynew)  # Построение кривой
    plt.grid(True)
    plt.title('Lagr')
    plt.savefig('lagr.png')
    plt.show()
