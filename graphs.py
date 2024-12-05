import numpy as np
import matplotlib.pyplot as plt


def test_tochnost(E1, E2, tol):
    """Проверка точности: если разница между E1 и E2 меньше точности, возвращает True."""
    return abs(E1 - E2) <= tol


def ekscentr_anom_iter(Ei, M, e, tol):
    """Рекурсивное вычисление эксцентрической аномалии с методом последовательных приближений."""
    EI = Ei
    Ei = M + e * np.sin(EI)

    if test_tochnost(Ei, EI, tol):
        return Ei
    else:
        return ekscentr_anom_iter(Ei, M, e, tol)


def true_anomaly(E, e):
    """Вычисление истинной аномалии с дополнительной проверкой значений."""
    try:
        return np.pi + 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan((np.pi + E) / 2))
    except:
        # Возврат NaN, если возникли проблемы с вычислением
        return np.nan


# Ввод параметров
rp = float(input("Введите радиус перигея (км): "))
ra = float(input("Введите радиус апогея (км): "))
R = float(input("Введите радиус планеты (км): "))
tol = float(input("Введите точность (например 0.0001): "))

# Вычисление эксцентриситета
e = (ra - rp) / (ra + rp + R * 2)
a = (ra + rp) / 2  # Большая полуось (среднее значение)

# Гравитационный параметр Земли (км^3/с^2)
mu = 398600

# Массивы для хранения данных
time_values = []
M_values = []
E_values = []
An_values = []
r_values = []  # Радиус
v_values = []  # Скорость

# Генерация данных
for t in np.linspace(0, 400, 401):
    M = 2 * np.pi * (t / 400)  # Средняя аномалия
    E = ekscentr_anom_iter(M, M, e, tol)  # Эксцентрическая аномалия
    An = true_anomaly(E, e)  # Истинная аномалия

    # Вычисление радиуса и скорости
    r = a * (1 - e ** 2) / (1 + e * np.cos(E))  # Радиус
    v = np.sqrt(mu * (2 / r - 1 / a))  # Скорость

    # Добавление данных в списки
    time_values.append(t)  # Время в секундах
    M_values.append(M / np.pi)  # Нормализованная средняя аномалия
    E_values.append(E / np.pi)  # Нормализованная эксцентрическая аномалия
    An_values.append(An / np.pi)  # Нормализованная истинная аномалия
    r_values.append(r)  # Радиус
    v_values.append(v)  # Скорость

# Сохранение данных в файлы
np.savetxt("T.txt", time_values)
np.savetxt("M.txt", M_values)
np.savetxt("E.txt", E_values)
np.savetxt("An.txt", An_values)
np.savetxt("r.txt", r_values)
np.savetxt("v.txt", v_values)

# Построение графиков
plt.figure(figsize=(12, 8))

# График истинной аномалии
plt.subplot(2, 2, 1)
plt.plot(time_values, An_values, label="Истинная аномалия (An)", color="blue")
plt.xlabel("Время (секунды)")
plt.ylabel("Истинная аномалия (An)")
plt.title("Зависимость истинной аномалии от времени")
plt.grid(True)

# График эксцентрической аномалии
plt.subplot(2, 2, 2)
plt.plot(time_values, E_values, label="Эксцентрическая аномалия (E)", color="green")
plt.xlabel("Время (секунды)")
plt.ylabel("Эксцентрическая аномалия (E)")
plt.title("Зависимость эксцентрической аномалии от времени")
plt.grid(True)

# График радиуса
plt.subplot(2, 2, 3)
plt.plot(time_values, r_values, label="Радиус (r)", color="purple")
plt.xlabel("Время (секунды)")
plt.ylabel("Радиус (км)")
plt.title("Зависимость радиуса от времени")
plt.legend()
plt.grid(True)

# График скорости
plt.subplot(2, 2, 4)
plt.plot(time_values, v_values, label="Скорость (v)", color="orange")
plt.xlabel("Время (секунды)")
plt.ylabel("Скорость (км/с)")
plt.title("Зависимость скорости от времени")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
