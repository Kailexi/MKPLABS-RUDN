import numpy as np
import matplotlib.pyplot as plt


def kepler_iteration_golden_section(M, e, a, b, tol, max_iter=1000):
    """
    Итерационный метод для нахождения корня уравнения Кеплера методом золотого сечения.

    :param M: Средняя аномалия (радианы)
    :param e: Эксцентриситет (0 <= e < 1)
    :param a: Левая граница интервала
    :param b: Правая граница интервала
    :param tol: Допустимая ошибка
    :param max_iter: Максимальное число итераций
    :return: Найденное значение эксцентрической аномалии E, список ошибок
    """
    # Функция для уравнения Кеплера
    def f(E):
        return E - e * np.sin(E) - M

    # Проверка начальных условий
    if f(a) * f(b) > 0:
        raise ValueError("Корень не лежит в заданном интервале [a, b]!")

    # Пропорция золотого сечения
    phi = (np.sqrt(5) + 1) / 2

    # Итерации
    errors = []
    for i in range(max_iter):
        # Определяем новые точки внутри интервала
        x1 = b - phi * (b - a)
        x2 = a + phi * (b - a)

        # Вычисление значений функции в этих точках
        f1, f2 = f(x1), f(x2)

        # Проверка точности
        if abs(b - a) < tol:
            return (a + b) / 2, errors

        # Сохраняем текущую ошибку
        errors.append(abs(f((a + b) / 2)))

        # Уточняем границы интервала
        if f1 * f2 < 0:
            # Корень находится между x1 и x2
            if f(a) * f1 < 0:
                b = x1
            else:
                a = x2
        else:
            if f1 == 0:
                return x1, errors
            if f2 == 0:
                return x2, errors
            # Корень находится между a и x1 или между x2 и b
            if f(a) * f1 < 0:
                b = x1
            else:
                a = x2

    raise RuntimeError("Метод не сошелся за заданное число итераций")



def kepler_iteration(M, e, a, b, x0, tol, max_iter=100):

    """
    Итерационный метод для нахождения корня уравнения Кеплера.

    :param M: Средняя аномалия (радианы)
    :param e: Эксцентриситет (0 <= e < 1)
    :param a: Левая граница интервала
    :param b: Правая граница интервала
    :param x0: Начальное приближение
    :param tol: Допустимая ошибка
    :param max_iter: Максимальное число итераций
    :return: Найденное значение эксцентрической аномалии E, список ошибок
    """


    # Функция для уравнения Кеплера
    def f(E):
        return E - e * np.sin(E) - M

    # Проверка начальных условий
    if f(a) * f(b) > 0:
        raise ValueError("Корень не лежит в заданном интервале [a, b]!")

    # Итерации
    x_prev = x0

    errors = []

    for i in range(max_iter):
        # Вычисление нового приближения
        x_next = (a + b) / 2

        # Сохранение ошибки
        errors.append(abs(f(x_next)))

        # Проверка точности
        if abs(f(x_next)) < tol:
            return x_next, errors

        # Обновление границ интервала
        if f(a) * f(x_next) < 0:

            b = x_next

        else:

            a = x_next


        x_prev = x_next  # Обновление начального приближения

    raise RuntimeError("Метод не сошелся за заданное число итераций")


M = float(input("Введите среднюю аномалию M (в радианах): "))
e = float(input("Введите эксцентриситет e (0 <= e < 1): "))
a = float(input("Введите левую границу интервала a: "))
b = float(input("Введите правую границу интервала b: "))
x0 = float(input("Введите начальное приближение x0: "))
tol = float(input("Введите точность ε: "))


try:

    E, errors = kepler_iteration(M, e, a, b, x0, tol)

    # Построение графика сходимости для метода половинного деления
    plt.figure(figsize=(8, 5))

    plt.plot(errors, marker="o", label="Ошибка на итерации при методе половинного деления")
    print(f"Найденное значение эксцентрической аномалии методом половинного деления E: {E:.6f}")


    E, errors = kepler_iteration_golden_section(M, e, a, b, tol)


    plt.plot(errors, marker='o', label="Ошибка на итерации при методе золотого сечения")

    plt.yscale('log')
    plt.xlabel("Итерация")
    plt.ylabel("Ошибка")

    plt.legend()

    print(f"Найденное значение эксцентрической аномалии методом золотого сечения E: {E:.6f}")

    # Построение графика сходимости для золотого сечения
    plt.title("Сравнения ошибкок на итерации для золотого сечения и метода половинного деления")
    plt.grid(True)
    plt.show()


except (ValueError, RuntimeError) as e:
    print(f"Ошибка: {e}")
