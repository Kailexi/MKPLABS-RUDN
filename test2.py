import numpy as np
import matplotlib.pyplot as plt


def test_tochnost(E1, E2, tochnost):
    """Проверка точности приближения эксцентрической аномалии."""
    return abs(E1 - E2) <= tochnost


def ekscentr_anom_iter(Ei, M, e, tochnost):
    """Итерационный метод для вычисления эксцентрической аномалии."""
    EI = Ei
    Ei = M + e * np.sin(EI)

    if test_tochnost(Ei, EI, tochnost):
        return Ei
    else:
        return ekscentr_anom_iter(Ei, M, e, tochnost)


# Вводные данные
rp = float(input("Введите радиус перигея: "))
ra = float(input("Введите радиус апогея: "))
R = float(input("Введите радиус планеты: "))
Mm = float(input("Введите массу планеты (например, 5.97): "))
stepen = int(input("Введите степень массы планеты (например, 24 для 10^24): "))
tochnost = float(input("Введите точность (например, 0.0001): "))

Mm = Mm * 10 ** stepen  # Приведение массы к полному значению


# Расчёт эксцентриситета
e = (ra - rp) / (ra + rp + R * 2)

# Константы
G = 6.67428e-20  # Гравитационная постоянная (км³/кг/с²)
p = 0.5 * (ra + rp + R * 2) * (1 - e ** 2)

# Временные параметры
time_steps = 400
time_values = np.linspace(0, 1, time_steps)

# Списки для хранения данных
M_values, E_values, An_values = [], [], []
r_values, Vr_values, Vn_values, V_values = [], [], [], []

# Основной цикл по времени
for t in time_values:
    M = 2 * np.pi * t  # Средняя аномалия
    E = ekscentr_anom_iter(M, M, e, tochnost)  # Эксцентрическая аномалия

    # Истинная аномалия через sin и cos
    sin_An = np.sqrt(1 - e ** 2) * np.sin(E) / (1 - e * np.cos(E))
    cos_An = (np.cos(E) - e) / (1 - e * np.cos(E))
    An = np.arctan2(sin_An, cos_An)  # atan2 учитывает квадранты

    # Приведение An к диапазону [0, 2π]
    if An < 0:
        An += 2 * np.pi

    r = p / (1 + e * np.cos(An))  # Радиус орбиты
    Vr = np.sqrt(G * Mm / p) * e * np.sin(An)  # Радиальная скорость
    Vn = np.sqrt(G * Mm / p) * (1 + e * np.cos(An))  # Тангенциальная скорость
    V = np.sqrt(Vn ** 2 + Vr ** 2)  # Полная скорость

    # Сохранение данных
    M_values.append(M / np.pi)
    E_values.append(E / np.pi)
    An_values.append(An / np.pi)
    r_values.append(r)
    Vr_values.append(Vr)
    Vn_values.append(Vn)
    V_values.append(V)

# Построение графиков

# График аномалий
plt.figure(figsize=(10, 6))
plt.plot(time_values, M_values, label='Средняя аномалия M(t)', color='b')
plt.plot(time_values, E_values, label='Эксцентрическая аномалия E(t)', color='g')
plt.plot(time_values, An_values, label='Истинная аномалия An(t)', color='r')
plt.xlabel('Время (t)')
plt.ylabel('Аномалия (в долях π)')
plt.legend()
plt.title('Зависимость аномалий от времени')
plt.grid(True)
plt.show()

# График радиуса
plt.figure(figsize=(10, 6))
plt.plot(time_values, r_values, label='Радиус r(t)', color='m')
plt.xlabel('Время (t)')
plt.ylabel('Радиус орбиты (км)')
plt.title('Зависимость радиуса от времени')
plt.grid(True)
plt.show()

# Графики скоростей
plt.figure(figsize=(10, 6))
plt.plot(time_values, Vr_values, label='Радиальная скорость Vr(t)', color='c')
plt.plot(time_values, Vn_values, label='Тангенциальная скорость Vn(t)', color='y')
plt.plot(time_values, V_values, label='Полная скорость V(t)', color='k')
plt.xlabel('Время (t)')
plt.ylabel('Скорость (км/с)')
plt.legend()
plt.title('Зависимость скоростей от времени')
plt.grid(True)
plt.show()
