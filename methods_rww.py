import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Gravitational constant, m^3/(kg·s^2)
M_earth = 5.972e24  # Earth's mass, kg
r_a_km = 46071  # Apogee, km
r_p_km = 6971  # Perigee, km

# Convert distances from kilometers to meters
r_a = r_a_km * 1e3
r_p = r_p_km * 1e3

# 1. Calculate semi-major axis and eccentricity
a = (r_a + r_p) / 2  # Semi-major axis, m
e = (r_a - r_p) / (r_a + r_p)  # Eccentricity

# 2. Calculate mean motion and orbital period
n = np.sqrt(G * M_earth / a**3)  # Mean motion, rad/s
T = 2 * np.pi / n  # Orbital period, s
T_days = T / (60 * 60 * 24)  # Orbital period in days

# Time intervals from 0 to T with a smaller step of 1 minute (60 seconds)
time_intervals = np.arange(0, T, 60)

# Mean anomaly M for each time step
M = n * time_intervals

# Kepler's equation solvers
def solve_kepler(M, e, tol=1e-6, max_iter=1000):
    E = M  # Initial guess
    for _ in range(max_iter):
        delta_E = (M - (E - e * np.sin(E))) / (1 - e * np.cos(E))
        E += delta_E
        if abs(delta_E) < tol:
            break
    return E

def iteration_method(M, e, tol=1e-6, max_iter=1000):
    E = M  # Initial guess
    for _ in range(max_iter):
        E_new = M + e * np.sin(E)
        if abs(E_new - E) < tol:
            return E_new
        E = E_new
    return E

def bisection_method(M, e, tol=1e-6, max_iter=1000):
    E_min, E_max = M - 2 * np.pi, M + 2 * np.pi
    for _ in range(max_iter):
        E_mid = (E_min + E_max) / 2
        f_mid = E_mid - e * np.sin(E_mid) - M
        if abs(f_mid) < tol:
            return E_mid
        if f_mid > 0:
            E_max = E_mid
        else:
            E_min = E_mid
    return (E_min + E_max) / 2

def golden_section_method(M, e, tol=1e-6, max_iter=1000):
    phi = (1 + np.sqrt(5)) / 2  # Golden ratio
    A = M - 2
    B = M + 2
    for _ in range(max_iter):
        C = A + (B - A) / phi
        f_C = C - e * np.sin(C) - M
        if abs(f_C) < tol:
            return C
        if (A - e * np.sin(A) - M) * f_C < 0:
            B = C
        else:
            A = C
    return (A + B) / 2

# Methods for solving Kepler's equation
methods = {
    "Метод Ньютона": solve_kepler,
    "Половинного Деления": bisection_method,
    "Золотого Сечения": golden_section_method,
    "Итераций": iteration_method,
}

# Plot anomalies for each method
for method_name, method in methods.items():
    E_values = np.array([method(Mi, e) for Mi in M])

    # True anomaly (v) from eccentric anomaly (E)
    v_values = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E_values / 2))
    v_values = np.where(v_values < 0, v_values + 2 * np.pi, v_values)

    # Plotting
    plt.figure(figsize=(12, 8))
    plt.plot(time_intervals / (60 * 60), M, label='Средняя Аномалия (M)', color='red', linewidth=1.5)
    plt.plot(time_intervals / (60 * 60), E_values, label='Эксцентрическая Аномалия (E)', color='green', linewidth=1.5)
    plt.plot(time_intervals / (60 * 60), v_values, label='Истинная Аномалия (v)', color='blue', linewidth=1.5)

    plt.title(f'Аномалии для метода:  {method_name}', fontsize=14)
    plt.xlabel('Время (часы)', fontsize=12)
    plt.ylabel('Аномалия (радианы)', fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.7)
    plt.legend()
    plt.tight_layout()
    plt.show()
