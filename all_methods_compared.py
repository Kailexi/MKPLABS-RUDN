
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
    """Метод половинного деления."""
    def f(E):
        return E - e * np.sin(E) - M

    errors = []
    for _ in range(max_iter):
        x_next = (a + b) / 2
        errors.append(abs(f(x_next)))

        if abs(f(x_next)) < tol:
            return x_next, errors

        if f(a) * f(x_next) < 0:
            b = x_next
        else:
            a = x_next

    raise RuntimeError("Метод не сошелся за заданное число итераций")


def kepler_newton(M, e, x0, tol, max_iter=100):
    """Метод Ньютона."""
    def f(E):
        return E - e * np.sin(E) - M

    def f_prime(E):
        return 1 - e * np.cos(E)

    errors = []
    x_prev = x0

    for _ in range(max_iter):
        x_next = x_prev - f(x_prev) / f_prime(x_prev)
        errors.append(abs(f(x_next)))

        if abs(f(x_next)) < tol:
            return x_next, errors

        x_prev = x_next

    raise RuntimeError("Метод не сошелся за заданное число итераций")


def kepler_fixed_point(M, e, x0, tol, max_iter=100):
    """Метод последовательных приближений."""
    def f(E):
        return M + e * np.sin(E)

    errors = []
    x_prev = x0

    for _ in range(max_iter):
        x_next = f(x_prev)
        errors.append(abs(x_next - x_prev))

        if abs(x_next - x_prev) < tol:
            return x_next, errors

        x_prev = x_next

    raise RuntimeError("Метод не сошелся за заданное число итераций")


def plot_solution(M, e, E, method_name):
    """Построение графика функции и вертикальной линии в точке корня."""
    E_values = np.linspace(a, b, 1000)
    f_values = E_values - e * np.sin(E_values) - M

    plt.plot(E_values, f_values, label="f(E)", color="blue")
    plt.axvline(x=E, color="red", linestyle="--", label=f"Найденное E = {E:.6f}")
    plt.axhline(y=0, color="black", linewidth=0.5)
    plt.title(f"Метод: {method_name}")
    plt.xlabel("E")
    plt.ylabel("f(E)")
    plt.grid(True)
    plt.legend()


def plot_errors(errors, method_name, ax):
    """Построение графика ошибки сходимости."""
    ax.plot(errors, marker="o", label=method_name)
    ax.set_yscale("log")
    ax.set_title(f"Ошибка: {method_name}")
    ax.set_xlabel("Итерация")
    ax.set_ylabel("Ошибка")
    ax.grid(True)


# Ввод параметров
M = float(input("Введите среднюю аномалию M (в радианах): "))
e = float(input("Введите эксцентриситет e (0 <= e < 1): "))
a = float(input("Введите левую границу интервала a: "))
b = float(input("Введите правую границу интервала b: "))
x0 = float(input("Введите начальное приближение x0: "))
tol = float(input("Введите точность ε: "))

try:
    # Метод половинного деления
    E_bisection, errors_bisection = kepler_iteration(M, e, a, b, x0, tol)
    print(f"Метод половинного деления: E = {E_bisection:.6f}")

    # Метод золотого сечения
    E_golden, errors_golden = kepler_iteration_golden_section(M, e, a, b, tol)
    print(f"Метод золотого сечения: E = {E_golden:.6f}")

    # Метод Ньютона
    E_newton, errors_newton = kepler_newton(M, e, x0, tol)
    print(f"Метод Ньютона: E = {E_newton:.6f}")

    # Метод последовательных приближений
    E_fixed_point, errors_fixed_point = kepler_fixed_point(M, e, x0, tol)
    print(f"Метод последовательных приближений: E = {E_fixed_point:.6f}")

    # Построение графиков решений
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Половинное деление
    plt.sca(axes[0, 0])
    plot_solution(M, e, E_bisection, "Метод половинного деления")

    # Золотое сечение
    plt.sca(axes[0, 1])
    plot_solution(M, e, E_golden, "Метод золотого сечения")

    # Метод Ньютона
    plt.sca(axes[1, 0])
    plot_solution(M, e, E_newton, "Метод Ньютона")

    # Последовательные приближения
    plt.sca(axes[1, 1])
    plot_solution(M, e, E_fixed_point, "Метод последовательных приближений")

    plt.tight_layout()
    plt.show()

    # Построение графиков ошибок
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Половинное деление
    plot_errors(errors_bisection, "Половинное деление", axes[0, 0])

    # Золотое сечение
    plot_errors(errors_golden, "Золотое сечение", axes[0, 1])

    # Метод Ньютона
    plot_errors(errors_newton, "Метод Ньютона", axes[1, 0])

    # Последовательные приближения
    plot_errors(errors_fixed_point, "Последовательные приближения", axes[1, 1])

    plt.tight_layout()
    plt.show()

except (ValueError, RuntimeError) as error:
    print(f"Ошибка: {error}")
