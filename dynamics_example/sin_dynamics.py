import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.interpolate import CubicSpline

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)

sys.path.insert(0, parent_dir + '/src/')

from propagation_solver import ModePropagator

# --- 1. ПАРАМЕТРЫ И КЭШИРОВАНИЕ ---
params = {
    'wavelength': 1.55, 'n_clad': 1.444, 'n_core': 3.4755,
    'W': 2.0, 'H': 0.22, 'd_xi': 0.05, 'd_eta': 0.05,
    'delta_u': 4.0, 'delta_d': 4.0, 'delta_l': 4.0, 'delta_r': 4.0 
}
NPML_l = [20, 20, 20, 20]

# Новое имя файла! Максимальная кривизна здесь ~0.16, кэш должен быть шире (до 0.2)
CACHE_FILE = "wg_cache_sinusoid_02.npz" 
propagator = ModePropagator(wg_params=params, num_modes=2, NPML=NPML_l, d_kappa=0.001)

if not os.path.exists(CACHE_FILE):
    k_array = np.linspace(-0.17, 0.17, 35)
    propagator.calculate_and_save_cache(k_array, filename=CACHE_FILE)

propagator.load_cache_and_build_splines(filename=CACHE_FILE)


# --- 2. ГЕОМЕТРИЯ СИНУСОИДАЛЬНОГО ВОЛНОВОДА ---
# Параметры из твоего документа
x0 = 100
A0 = 10

# Создаем плотную сетку по x
x_arr = np.linspace(0, x0, 1000)

# Формулы производных из документа (для уравнения 33 и 34)
y_prime = A0 * (4 * np.pi / x0) * np.cos(4 * np.pi * x_arr / x0)
y_prime2 = -A0 * (4 * np.pi / x0)**2 * np.sin(4 * np.pi * x_arr / x0)

# Вычисляем координату s(x) через интеграл 
s_integrand = np.sqrt(1 + y_prime**2)
s_arr = integrate.cumulative_trapezoid(s_integrand, x_arr, initial=0.0)
L = s_arr[-1] # Полная длина волновода L = s(x0)

# Вычисляем кривизну kappa(x) 
kappa_arr = y_prime2 / (1 + y_prime**2)**1.5

# МАГИЯ: Строим кубический сплайн kappa(s) по нашим массивам!
# Сплайн сам умеет брать точные производные первого и второго порядка
kappa_spline = CubicSpline(s_arr, kappa_arr)

def k_func(s):
    return kappa_spline(s)

def dk_func(s):
    return kappa_spline(s, 1) # Первая производная по s

def d2k_func(s):
    return kappa_spline(s, 2) # Вторая производная по s


# --- 3. РЕШЕНИЕ ДИНАМИКИ ---
c0 = np.array([1.0, 0.0], dtype=complex)
sol = propagator.solve_dynamics([0, L], k_func, dk_func, d2k_func, c0)


# --- 4. ОТРИСОВКА ---
s_vals = np.linspace(0, L, 500)
Y_vals = sol.sol(s_vals)

P1 = np.abs(Y_vals[0, :])**2
P2 = np.abs(Y_vals[1, :])**2
P_total = P1 + P2

# Восстанавливаем x и y для красивого графика
phi_vals = integrate.cumulative_trapezoid(k_func(s_vals), s_vals, initial=0.0)
x_sim = integrate.cumulative_trapezoid(np.cos(phi_vals), s_vals, initial=0.0)
y_sim = integrate.cumulative_trapezoid(np.sin(phi_vals), s_vals, initial=0.0)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))

# График 1: Физическая форма волновода (как на Рис. 5)
ax1.plot(x_sim, y_sim, 'g-', linewidth=4, label='Waveguide Core')
ax1.set_title('Physical Shape of the Sinusoidal Waveguide', fontsize=14)
ax1.set_xlabel(r'$x, \mu m$', fontsize=12)
ax1.set_ylabel(r'$y, \mu m$', fontsize=12)
ax1.set_aspect('equal')
ax1.grid(True)
ax1.legend()

# График 2: Функция кривизны (как на Рис. 6)
ax2.plot(s_vals, k_func(s_vals), 'b-', linewidth=3)
ax2.set_title(r'Curvature Function $\kappa(s)$', fontsize=14)
ax2.set_xlabel(r'$s, \mu m$', fontsize=12)
ax2.set_ylabel(r'$\kappa, \mu m^{-1}$', fontsize=12)
ax2.grid(True)

# График 3: Динамика перекачки энергии
ax3.plot(s_vals, P1, label='Mode 1 Power', color='blue', linewidth=2)
ax3.plot(s_vals, P2, label='Mode 2 Power', color='red', linewidth=2)
ax3.plot(s_vals, P_total, 'k--', label='Total Power', alpha=0.7)
ax3.set_title('Dynamics of Mode Coupling', fontsize=14)
ax3.set_xlabel(r'Propagation distance $s$, $\mu m$', fontsize=12)
ax3.set_ylabel('Normalized Power', fontsize=12)
ax3.set_ylim([0, 1.05])
ax3.legend(fontsize=12)
ax3.grid(True)

plt.tight_layout()
plt.show()