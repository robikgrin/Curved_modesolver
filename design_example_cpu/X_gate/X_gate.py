import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.integrate as integrate
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
from scipy.linalg import dft
from tqdm import tqdm

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)

sys.path.insert(0,  os.path.dirname(parent_dir) + '/src/')

from propagation_solver import ModePropagator

# =====================================================================
# 1. ИНИЦИАЛИЗАЦИЯ И ЗАГРУЗКА КЭША
# =====================================================================
params = {
    'wavelength': 1.55, 'n_clad': 1.444, 'n_core': 3.4755,
    'W': 1, 'H': 0.22, 'd_xi': 0.05, 'd_eta': 0.05,
    'delta_u': 4.0, 'delta_d': 4.0, 'delta_l': 4.0, 'delta_r': 4.0 
}

# Микрометровый масштаб
N_modes = 2
propagator = ModePropagator(wg_params=params, num_modes=N_modes, NPML=[20, 20, 20, 20], d_kappa=0.001)

CACHE_FILE = "wg_cache_micro.npz"
if not os.path.exists(CACHE_FILE):
    print("Создание нового кэша матриц (СИ - микрометры)...")
    k_array = np.linspace(-0.2, 0.2, 35) # Безопасный диапазон без сингулярностей
    propagator.calculate_and_save_cache(k_array, filename=CACHE_FILE)

print("Загрузка Мега-сплайна...")
propagator.load_cache_and_build_splines(filename=CACHE_FILE)


# =====================================================================
# 2. ПАРАМЕТРИЗАЦИЯ ВОЛНОВОДА (Даем больше свободы)
# =====================================================================
L = 100      # Заданная длина волновода
N_points = 30  # Число оптимизировочных точек (контрольных точек сплайна)

def get_curvature_funcs_spline(params, L):
    s_points = np.linspace(L/(N_points+1), L - L/(N_points+1), N_points)
    s_full = np.concatenate(([0.0], s_points, [L]))
    k_full = np.concatenate(([0.0], params, [0.0]))
    spline = CubicSpline(s_full, k_full, bc_type='clamped')
    return spline, spline.derivative(1), spline.derivative(2)

# ЦЕЛЕВАЯ ФУНКЦИЯ
U_target = np.array([[0, 1], [1, 0]], dtype=complex)

def objective_function(params):
    # Жесткий барьер для защиты от выхода за пределы кэша
    if np.max(np.abs(params)) > 0.15:
        return 10.0 + np.sum(params**2)

    k_func, dk_func, d2k_func = get_curvature_funcs_spline(params, L)
    
    M_sim = np.zeros((N_modes, N_modes), dtype=complex)
    
    # Прогоняем диффур N раз для каждого порта входа
    for i in range(N_modes):
        c0 = np.zeros(N_modes, dtype=complex)
        c0[i] = 1.0 
        
        sol = propagator.solve_dynamics([0, L], k_func, dk_func, d2k_func, c0, fast_mode=False)
        M_sim[:, i] = sol.y[0:N_modes, -1]
    
    # ==========================================================
    # 1. ОШИБКА КВАНТОВОГО ВЕНТИЛЯ (Fidelity Error)
    # ==========================================================
    trace_val = np.trace(np.conjugate(U_target.T) @ M_sim)
    fidelity = (np.abs(trace_val)**2) / (N_modes**2)
    error_fidelity = np.abs(1.0 - fidelity)
    
    # ==========================================================
    # 2. ФИЗИЧЕСКИЕ ПОТЕРИ (Insertion Loss Penalty)
    # ==========================================================
    # Считаем среднее пропускание устройства по всем портам
    transmission = np.sum(np.abs(M_sim)**2) / N_modes
    # Штраф за потери: 0.0 если 100% фотонов долетели, стремится к 1.0 при высоких потерях
    error_loss = np.abs(1.0 - transmission) 
    
    # ==========================================================
    # 3. ШТРАФ ЗА "ИЗВИЛИСТОСТЬ" (Fabrication Penalty)
    # ==========================================================
    penalty_smoothness = 0.001 * np.sum(params**2)
    
    # --- ВЕСОВЫЕ КОЭФФИЦИЕНТЫ ---
    # Мы можем управлять приоритетами алгоритма.
    # Если alpha_loss = 2.0, алгоритм будет яростно бороться за каждый фотон, 
    # даже немного жертвуя точностью фазы.
    alpha_loss = 0.1
    
    return error_fidelity + (alpha_loss * error_loss) + penalty_smoothness

# =====================================================================
# 4. АЛГОРИТМ ОБРАТНОГО ПРОЕКТИРОВАНИЯ
# =====================================================================
initial_guess = 0.05 * np.sin(np.linspace(0, 4 * np.pi, N_points))

max_iter = 1000
pbar = tqdm(total=max_iter, desc="Оптимизация формы (Nelder-Mead)", unit="итер")

def callback(x):
    pbar.update(1)

print("\nСтарт оптимизации")
res = minimize(
    objective_function, 
    initial_guess, 
    method='Nelder-Mead', 
    callback=callback,
    options={'maxiter': max_iter, 'xatol': 1e-4, 'fatol': 1e-4}
)
pbar.close()

print(f"\nОПТИМИЗАЦИЯ УСПЕШНА!")
print(f"Найденные параметры сплайна:\n{np.round(res.x, 4)}")


# =====================================================================
# ФИНАЛЬНЫЙ РАСЧЕТ И ОТРИСОВКА ИДЕАЛЬНОГО ВОЛНОВОДА
# =====================================================================
k_func, dk_func, d2k_func = get_curvature_funcs_spline(res.x, L)

# 1. Расчет полной матрицы передачи M_sim (для 3D гистограммы)
M_sim = np.zeros((N_modes, N_modes), dtype=complex)
for i in range(N_modes):
    c0_bench = np.zeros(N_modes, dtype=complex)
    c0_bench[i] = 1.0
    sol_m = propagator.solve_dynamics([0, L], k_func, dk_func, d2k_func, c0_bench, fast_mode=False)
    M_sim[:, i] = sol_m.y[0:N_modes, -1]

# 2. Детальный прогон динамики для анализа потерь (возбуждаем 1-ю моду)
c0_init = np.zeros(N_modes, dtype=complex)
c0_init[0] = 1.0
sol_final = propagator.solve_dynamics([0, L], k_func, dk_func, d2k_func, c0_init, fast_mode=False)

s_vals = np.linspace(0, L, 500)
Y_all = sol_final.sol(s_vals)

# Мощности каждой моды и суммарная физическая мощность
P_modes = np.abs(Y_all[0:N_modes, :])**2
P_total = np.sum(P_modes, axis=0)

trace_val_final = np.trace(np.conjugate(U_target.T) @ M_sim)
final_fidelity = (np.abs(trace_val_final)**2) / (N_modes**2)

final_transmission = np.sum(np.abs(M_sim)**2) / N_modes

print(f"\n" + "="*40)
print(f" ФИНАЛЬНЫЙ ОТЧЕТ (N={N_modes} мод, L={L} мкм)")
print(f"========================================")
print(f"Квантовая Фиделити (Fidelity): {final_fidelity*100:.2f}%")
print(f"Среднее пропускание (Transmission): {final_transmission*100:.2f}%")
print(f"Суммарные потери (Insertion Loss): {-10 * np.log10(final_transmission):.3f} dB")
print(f"="*40 + "\n")

# --- МАТЕМАТИКА ПОТЕРЬ В dB ---
loss_db = -10 * np.log10(np.clip(P_total, 1e-12, 1.0))
ds = s_vals[1] - s_vals[0]
loss_rate = np.gradient(loss_db, ds)

# Теоретические радиационные потери (из мнимой части бета фундаментальной моды)
vec_data = propagator.mega_spline(k_func(s_vals))
beta_imag = np.abs(np.imag(vec_data[:, 0])) 
theoretical_rate = 2 * beta_imag * 4.343 # Коэффициент пересчета в dB/um

# Восстанавливаем физическую геометрию x(s), y(s)
phi_vals = integrate.cumulative_trapezoid(k_func(s_vals), s_vals, initial=0.0)
x_sim = integrate.cumulative_trapezoid(np.cos(phi_vals), s_vals, initial=0.0)
y_sim = integrate.cumulative_trapezoid(np.sin(phi_vals), s_vals, initial=0.0)

# =====================================================================
# ВИЗУАЛИЗАЦИЯ (Презентационный стиль)
# =====================================================================
plt.style.use('seaborn-v0_8-whitegrid') # Делаем графики чище
colors = plt.cm.viridis(np.linspace(0, 0.9, N_modes))

# 1. ГРАФИК 1: ФОРМА И ДИНАМИКА (Склеенный)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), dpi=200)

# Геометрия
ax1.plot(x_sim, y_sim, color='teal', linewidth=4, solid_capstyle='round')
ax1.set_title('Inverse Designed X gate Interferometer', fontsize=16, pad=15)
ax1.set_xlabel(r'$x, \mu m$', fontsize=12)
ax1.set_ylabel(r'$y, \mu m$', fontsize=12)
ax1.set_aspect('equal')
ax1.grid(True, linestyle='--', alpha=0.7)

# Динамика мощностей
for i in range(N_modes):
    ax2.plot(s_vals, P_modes[i], color=colors[i], linewidth=2.5, label=f'Mode {i+1}')
ax2.plot(s_vals, P_total, color='black', linestyle='--', linewidth=2, alpha=0.8, label='Total Power')
ax2.set_title(f'Power Dynamics (Exciting Mode 1)', fontsize=14, pad=10)
ax2.set_xlabel(r'Propagation distance $s$, $\mu m$', fontsize=12)
ax2.set_ylabel('Normalized Power', fontsize=12)
ax2.set_ylim([-0.05, 1.05])
ax2.legend(loc='upper right', frameon=True, fontsize=10)
ax2.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig("geometry_and_dynamics.png", dpi=400, bbox_inches='tight')
plt.show()

# 2. ГРАФИК 2: ГЛУБОКИЙ АНАЛИЗ ПОТЕРЬ
fig_loss, (ax_l1, ax_l2) = plt.subplots(2, 1, figsize=(10, 8), dpi=200, sharex=True)

# Накопленные потери
ax_l1.plot(s_vals, loss_db, color='darkred', linewidth=3)
ax_l1.set_ylabel('Cumulative Loss [dB]', fontsize=12)
ax_l1.set_title(f'Insertion Loss Profile (Total IL: {loss_db[-1]:.3f} dB)', fontsize=14)
ax_l1.grid(True, linestyle='--', alpha=0.7)

# Локальная интенсивность (Hotspots)
ax_l2.fill_between(s_vals, loss_rate, color='coral', alpha=0.4, label='Total Simulated Loss Rate')
ax_l2.plot(s_vals, loss_rate, color='orangered', linewidth=2)
ax_l2.plot(s_vals, theoretical_rate, color='darkred', linestyle='--', linewidth=2, alpha=0.8, label='Pure Bending Loss (PML theory)')

ax_l2.set_ylabel(r'Loss Rate [dB/$\mu m$]', fontsize=12)
ax_l2.set_xlabel(r'Propagation distance $s$, $\mu m$', fontsize=12)
ax_l2.set_title('Loss Hotspots (Transition vs Bending)', fontsize=14)
ax_l2.legend(loc='upper right', frameon=True)
ax_l2.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig("loss_analysis.png", dpi=400, bbox_inches='tight')
plt.show()

# 3. ГРАФИК 3: СРАВНЕНИЕ МАТРИЦ (3D)
def plot_complex_matrix_3d(ax, matrix, title):
    N = matrix.shape[0]
    x, y = np.meshgrid(np.arange(N), np.arange(N))
    x, y = x.flatten(), y.flatten()
    z = np.zeros_like(x)
    
    magnitudes = np.abs(matrix).flatten()
    phases = np.angle(matrix).flatten()
    
    norm_phases = (phases + np.pi) / (2 * np.pi)
    cols = cm.hsv(norm_phases)
    
    dx = dy = 0.6
    ax.bar3d(x - dx/2, y - dy/2, z, dx, dy, magnitudes, color=cols, shade=True, alpha=0.9)
    
    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel('Input Mode', labelpad=10)
    ax.set_ylabel('Output Mode', labelpad=10)
    ax.set_zlabel(r'Amplitude $|U_{ij}|$', labelpad=10)
    ax.set_xticks(np.arange(N))
    ax.set_yticks(np.arange(N))
    ax.view_init(elev=30, azim=45) # Делаем угол обзора чуть красивее

fig_3d = plt.figure(figsize=(14, 6), dpi=200)

ax1_3d = fig_3d.add_subplot(121, projection='3d')
plot_complex_matrix_3d(ax1_3d, U_target, r'Target Unitary $U_{target}$ (DFT)')

ax2_3d = fig_3d.add_subplot(122, projection='3d')
plot_complex_matrix_3d(ax2_3d, M_sim, r'Experimental Unitary $M_{sim}$')

sm = cm.ScalarMappable(cmap=cm.hsv, norm=plt.Normalize(vmin=-np.pi, vmax=np.pi))
cbar_ax = fig_3d.add_axes([0.92, 0.25, 0.015, 0.5]) # Тонкий и аккуратный colorbar
cbar = fig_3d.colorbar(sm, cax=cbar_ax)
cbar.set_label(r'Phase $\arg(U_{ij})$', fontsize=12)
cbar.set_ticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
cbar.set_ticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])

plt.subplots_adjust(wspace=0.1)
plt.savefig('unitary_matrices_3d.png', dpi=400, bbox_inches='tight')
plt.show()


# =====================================================================
# 4. ГРАФИК 4: ФУНКЦИЯ КРИВИЗНЫ (Curvature Profile)
# =====================================================================
fig_curv, ax_c = plt.subplots(figsize=(10, 5), dpi=200)

# Рисуем саму функцию кривизны
kappa_opt = k_func(s_vals)
ax_c.plot(s_vals, kappa_opt, color='purple', linewidth=2.5, label=r'Optimized $\kappa(s)$')

# Добавляем заливку для наглядности (показывает изгибы влево/вправо)
ax_c.fill_between(s_vals, kappa_opt, 0, where=(kappa_opt > 0), color='purple', alpha=0.2)
ax_c.fill_between(s_vals, kappa_opt, 0, where=(kappa_opt <= 0), color='indigo', alpha=0.2)

# Линия нулевой кривизны (прямой волновод)
ax_c.axhline(0, color='black', linewidth=1, linestyle='--')

ax_c.set_title('Optimized Curvature Profile', fontsize=15, pad=10)
ax_c.set_xlabel(r'Propagation distance $s$, $\mu m$', fontsize=12)
ax_c.set_ylabel(r'Curvature $\kappa(s)$, $\mu m^{-1}$', fontsize=12)
ax_c.grid(True, linestyle='--', alpha=0.7)
ax_c.legend(loc='upper right', frameon=True, fontsize=11)

plt.tight_layout()
plt.savefig("optimized_curvature.png", dpi=400, bbox_inches='tight')
plt.show()