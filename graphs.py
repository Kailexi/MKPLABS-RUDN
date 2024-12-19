import numpy as np
import matplotlib.pyplot as plt
import csv

# Константы
G = 6.67430e-11  # гравитационная постоянная (м^3/кг/с^2)

# Ввод параметров
M_Earth = float(input("Введите массу Земли (кг): "))  # масса Земли (кг)
r_a_km = float(input("Введите радиус апоцентра (км): "))  # апоцентр (км)
r_p_km = float(input("Введите радиус перицентра (км): "))  # перицентр (км)

r_a = r_a_km * 1e3  # апоцентр в метрах
r_p = r_p_km * 1e3  # перицентр в метрах

# Полуось и эксцентриситет
a = (r_a + r_p) / 2  # большая полуось
e = (r_a - r_p) / (r_a + r_p)  # эксцентриситет

# Среднее движение и период
mu = G * M_Earth
n = np.sqrt(mu / a ** 3)  # среднее движение
T = 2 * np.pi / n  # орбитальный период
time_intervals = np.linspace(0, T, 1000)  # временные интервалы


# Метод Ньютона для решения уравнения Кеплера
def newton_method(M, e, epsilon=1e-6, max_it=1000):
    E = M
    for _ in range(max_it):
        f = E - e * np.sin(E) - M
        df = 1 - e * np.cos(E)
        dE = -f / df
        E += dE
        if abs(dE) < epsilon:
            break
    return E


# Расчет параметров орбиты
def calculate_orbit_params(time_intervals, e, a, n, mu, epsilon=1e-6, max_it=1000):
    p = a * (1 - e ** 2)  # фокальный параметр
    r_values, Vr_values, Vt_values, V_values = [], [], [], []
    for t in time_intervals:
        M = n * t
        E = newton_method(M, e, epsilon, max_it)
        v = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
        r = p / (1 + e * np.cos(v))
        Vr = np.sqrt(mu / p) * e * np.sin(v)
        Vt = np.sqrt(mu / p) * (1 + e * np.cos(v))
        V = np.sqrt(Vr ** 2 + Vt ** 2)
        r_values.append(r)
        Vr_values.append(Vr)
        Vt_values.append(Vt)
        V_values.append(V)
    return np.array(r_values), np.array(Vr_values), np.array(Vt_values), np.array(V_values)


# Выполнение расчетов
r, Vr, Vt, V = calculate_orbit_params(time_intervals, e, a, n, mu)

# Найдем максимальное и минимальное значение радиальной скорости
Vr_max = np.max(Vr)
Vr_min = np.min(Vr)
V_at_Vr_max = V[np.argmax(Vr)]
V_at_Vr_min = V[np.argmin(Vr)]

# Найдем максимальное и минимальное значение трансверсальной скорости
Vt_max = np.max(Vt)
Vt_min = np.min(Vt)
V_at_Vt_max = V[np.argmax(Vt)]
V_at_Vt_min = V[np.argmin(Vt)]

# Найдем максимальный и минимальный радиус
R_max = np.max(r)
R_min = np.min(r)

# Результаты
print(f"Максимальный радиус: {R_max / 1e3:.2f} км")
print(f"Минимальный радиус: {R_min / 1e3:.2f} км")
print(f"Максимальная радиальная скорость: {Vr_max / 1e3:.2f} км/с")
print(f"Полная скорость при максимальной радиальной: {V_at_Vr_max / 1e3:.2f} км/с")
print(f"Минимальная радиальная скорость: {Vr_min / 1e3:.2f} км/с")
print(f"Полная скорость при минимальной радиальной: {V_at_Vr_min / 1e3:.2f} км/с")
print(f"Максимальная трансверсальная скорость: {Vt_max / 1e3:.2f} км/с")
print(f"Полная скорость при максимальной трансверсальной: {V_at_Vt_max / 1e3:.2f} км/с")
print(f"Минимальная трансверсальная скорость: {Vt_min / 1e3:.2f} км/с")
print(f"Полная скорость при минимальной трансверсальной: {V_at_Vt_min / 1e3:.2f} км/с")

# Запись данных в CSV файл
with open('orbit_data.csv', 'w', newline='') as csvfile:
    fieldnames = ['Time (hours)', 'Radius (m)', 'Vr (m/s)', 'Vt (m/s)', 'V (m/s)']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for t, r_val, Vr_val, Vt_val, V_val in zip(time_intervals / 3600, r, Vr, Vt, V):
        writer.writerow({
            'Time (hours)': t,
            'Radius (m)': r_val,
            'Vr (m/s)': Vr_val,
            'Vt (m/s)': Vt_val,
            'V (m/s)': V_val
        })

print("Данные успешно записаны в 'orbit_data.csv'.")

# Построение графиков
plt.figure(figsize=(10, 6))
plt.plot(time_intervals / 3600, r, label='Радиус-вектор (м)')
plt.title('Радиус-вектор объекта')
plt.xlabel('Время (часы)')
plt.ylabel('Радиус-вектор (м)')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(time_intervals / 3600, Vr / 1e3, label='Радиальная скорость (км/с)')
plt.plot(time_intervals / 3600, Vt / 1e3, label='Трансверсальная скорость (км/с)')
plt.plot(time_intervals / 3600, V / 1e3, label='Полная скорость (км/с)')
plt.title('Скорости объекта')
plt.xlabel('Время (часы)')
plt.ylabel('Скорость (км/с)')
plt.grid(True)
plt.legend()
plt.show()
