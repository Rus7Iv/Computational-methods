import numpy

from matplotlib import pyplot
from scipy import interpolate

m = int(input("Введите количество букв в вашей фамилии (Метелев = 7): "))
k = int(input("Введите количество букв в вашем имени (Виталий = 7): "))

N = 3 * k + 2

x = numpy.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
y = numpy.array([0.2 * N, 0.3 * m, 0.5 * k, 0.6 * N,
                 0.7 * m, k, 0.8 * N, 1.2 * k, 1.3 * m, N])

# Метод Лагранжа
def lagrange(x, y, t):
    result = 0
    for j in range(len(y)):
        point1 = 1
        point2 = 1
        for i in range(len(x)):
            if i != j:
                point1 = point1 * (t - x[i])
                point2 = point2 * (x[j] - x[i])
        result = result + y[j] * point1 / point2
    return result


# Наилучшее среднеквадратичное приближение
def rootMeanSquareApproximation(x, y):
    currentX = numpy.linspace(min(x), max(x), num=numpy.size(x) * 100)
    coefficients = numpy.polyfit(x, y, 2)
    currentLine = numpy.polyval(coefficients, currentX)
    pyplot.title("Наилучшее среднеквадратичное приближение")
    pyplot.scatter(x, y)
    pyplot.scatter(currentX, currentLine, c='b', s=1)
    pyplot.xlim(min(x) - 0.05, max(x) + 0.05)
    pyplot.xticks(rotation=90)
    pyplot.tight_layout()
    pyplot.show()


def main():
    methodNumber = int(input("\n" +
                             "  1: МНОГОЧЛЕН ЛАГРАНЖА\n" +
                             "  2: ПАРАБОЛИЧЕСКИЙ СПЛАЙН\n" +
                             "  3: НАИЛУЧШЕЕ СРЕДНЕКВАДРАТИЧНОЕ ПРИБЛИЖЕНИЕ\n\n" +
                             " Выберите то, что нужно найти и построить: "))
    if methodNumber == 1:
        interpolatedX = numpy.linspace(numpy.min(x), numpy.max(x), 100)
        interpolatedY = [lagrange(x, y, i) for i in interpolatedX]
        pyplot.plot(interpolatedX, interpolatedY, 'r')
        pyplot.title("Многочлен Лагранжа")
        pyplot.grid(True)
        pyplot.show()

    if methodNumber == 2:
        tck = interpolate.splrep(x, y, s=0)
        ynew = interpolate.splev(x, tck, der=0)
        pyplot.title("Параболический сплайн")
        pyplot.plot(x, ynew, '-')
        pyplot.grid(True)
        pyplot.show()

    if methodNumber == 3:
        rootMeanSquareApproximation(x, y)


main()
