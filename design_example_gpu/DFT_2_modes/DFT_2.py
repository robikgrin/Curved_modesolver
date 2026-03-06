import os
import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import scipy.integrate as integrate
from scipy.linalg import dft
import time
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)

sys.path.insert(0, os.path.dirname(parent_dir) + '/gpu/')

from gpu_solver import TorchModePropagator

sys.path.insert(0, os.path.dirname(parent_dir) + '/src/')

from propagation_solver import ModePropagator

# =====================================================================
# НАСТРОЙКИ И ИНИЦИАЛИЗАЦИЯ
# =====================================================================
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Вычисления запущены на: {device}")

params = {
    'wavelength': 1.55, 'n_clad': 1.444, 'n_core': 3.4755,
    'W': 1, 'H': 0.22, 'd_xi': 0.05, 'd_eta': 0.05,
    'delta_u': 4.0, 'delta_d': 4.0, 'delta_l': 4.0, 'delta_r': 4.0 
}

CACHE_FILE = "wg_W1_N2.npz" 
L = 100.0        # Увеличим длину, чтобы дать свету место для перекачки 50/50
num_steps = 1000 # ОЧЕНЬ ВАЖНО: ds = 0.1 мкм. Только так RK4 выживет!
N_modes = 2

NPML_l = [20, 20, 20, 20]

propagator =  ModePropagator(wg_params=params, num_modes=N_modes, NPML=NPML_l, d_kappa=0.001)

if not os.path.exists(CACHE_FILE):
    k_array = np.linspace(-0.17, 0.17, 35)
    propagator.calculate_and_save_cache(k_array, filename=CACHE_FILE)

solver = TorchModePropagator(CACHE_FILE, device=device)

# --- ПАРАМЕТРИЗАЦИЯ ВОЛНОВОДА ---
N_ctrl = 20 # 60 точек достаточно для гибкости, 100 - это перебор

# НАЧИНАЕМ С ПЛАВНОЙ КРИВОЙ, А НЕ СО СЛУЧАЙНОГО ШУМА
initial_k = 0.05 * np.sin(np.linspace(0, 4 * np.pi, N_ctrl))
k_raw = nn.Parameter(torch.tensor(initial_k, dtype=torch.float64, device=device))

# Теперь оптимизируем И кривизну, И выходные фазы
optimizer = torch.optim.Adam([k_raw], lr=0.01)

# --- ЦЕЛЕВАЯ МАТРИЦА (DFT) ---
U_target_np = dft(N_modes) / np.sqrt(N_modes)
Target_M = torch.tensor(U_target_np, dtype=torch.complex128, device=device)

# =====================================================================
# ЦИКЛ ОБУЧЕНИЯ (Adam)
# =====================================================================
print("\nЗапуск Градиентного Обратного Проектирования (Adam)...")
start_time = time.time()
epochs = 800
loss_history, fidelity_history = [], []

for epoch in range(epochs):
    optimizer.zero_grad()
    
    # Ограничение кривизны
    k_ctrl = 0.19 * torch.tanh(k_raw)
    
    # Вызов RK4 PyTorch солвераß
    M_sim = solver(L, k_ctrl)
    
    # Считаем Фиделити НАПРЯМУЮ, без всяких дополнительных фаз
    trace_val = torch.trace(torch.matmul(Target_M.mH, M_sim))
    fidelity = (torch.abs(trace_val)**2) / (N_modes**2)
    error_fidelity = torch.abs(1.0 - fidelity)
    
    transmission = torch.sum(torch.abs(M_sim)**2) / N_modes
    error_loss = torch.abs(1.0 - transmission)
    
    smoothness_penalty = 0.00001 * torch.sum((k_ctrl[1:] - k_ctrl[:-1])**2)
    
    loss = error_fidelity + (0.5 * error_loss) + smoothness_penalty
    loss.backward()
    optimizer.step()
    
    loss_history.append(loss.item())
    fidelity_history.append(fidelity.item())
    
    if (epoch + 1) % 50 == 0:
        print(f"Итерация {epoch+1}/{epochs} | Loss: {loss.item():.4f} | Fidelity: {fidelity.item()*100:.2f}% | Transm: {transmission.item()*100:.2f}%")
print(f"Оптимизация завершена за {time.time() - start_time:.2f} секунд!")

# --- ФИНАЛЬНЫЙ АНАЛИЗ ---
with torch.no_grad():
    solver.eval() 
    k_ctrl_final = 0.2 * torch.tanh(k_raw)
    
    # Солвер теперь возвращает 2 значения, как мы писали в RK4:
    M_sim, k_dense = solver(L, k_ctrl_final)
    
    phase_diag = torch.exp(1j * out_phases).unsqueeze(1)
    M_corrected = phase_diag * M_sim
    
    trace_val = torch.trace(torch.matmul(Target_M.mH, M_corrected))
    final_fidelity = ((torch.abs(trace_val)**2) / (N_modes**2)).item()
    
    M_final = M_corrected.cpu().numpy()
    final_transmission = np.sum(np.abs(M_final)**2) / N_modes
    kappa_final = k_dense.cpu().numpy()

print(f"\n" + "="*40)
print(f" ФИНАЛЬНЫЙ ОТЧЕТ (PyTorch GPU + Nelder-Mead)")
print(f"========================================")
print(f"Квантовая Фиделити: {final_fidelity*100:.2f}%")
print(f"Среднее пропускание: {final_transmission*100:.2f}%")
print(f"Суммарные потери (IL): {-10 * np.log10(final_transmission):.3f} dB")
print(f"="*40 + "\n")
# --- ВОССТАНОВЛЕНИЕ ДИНАМИКИ НА CPU ДЛЯ ГРАФИКОВ ---

# Создаем сплайн из тензора kappa_final
spline_k = CubicSpline(s_vals, kappa_final)
dk_func = spline_k.derivative(1)
d2k_func = spline_k.derivative(2)

cpu_propagator = ModePropagator(wg_params=solver.N, num_modes=N_modes, NPML=[20, 20, 20, 20], d_kappa=0.001)
cpu_propagator.load_cache_and_build_splines(filename=CACHE_FILE)

c0_init = np.zeros(N_modes, dtype=complex)
c0_init[0] = 1.0 # Анализируем возбуждение первой моды
sol_final = cpu_propagator.solve_dynamics([0, L], spline_k, dk_func, d2k_func, c0_init, fast_mode=False)

Y_all = sol_final.sol(s_vals)
P_modes = np.abs(Y_all[0:N_modes, :])**2
P_total = np.sum(P_modes, axis=0)

# Математика потерь
loss_db = -10 * np.log10(np.clip(P_total, 1e-12, 1.0))
ds = s_vals[1] - s_vals[0]
loss_rate = np.gradient(loss_db, ds)

vec_data = cpu_propagator.mega_spline(kappa_final)
beta_imag = np.abs(np.imag(vec_data[:, 0])) 
theoretical_rate = 2 * beta_imag * 4.343

phi_vals = integrate.cumulative_trapezoid(kappa_final, s_vals, initial=0.0)
x_sim = integrate.cumulative_trapezoid(np.cos(phi_vals), s_vals, initial=0.0)
y_sim = integrate.cumulative_trapezoid(np.sin(phi_vals), s_vals, initial=0.0)

# =====================================================================
# ВИЗУАЛИЗАЦИЯ (Презентационный стиль)
# =====================================================================
plt.style.use('seaborn-v0_8-whitegrid')
colors = plt.cm.viridis(np.linspace(0, 0.9, N_modes))

# 1. Форма и Динамика
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), dpi=200)
ax1.plot(x_sim, y_sim, color='teal', linewidth=4, solid_capstyle='round')
ax1.set_title('Inverse Designed DFT Interferometer (PyTorch)', fontsize=16, pad=15)
ax1.set_xlabel(r'$x, \mu m$', fontsize=12); ax1.set_ylabel(r'$y, \mu m$', fontsize=12)
ax1.set_aspect('equal'); ax1.grid(True, linestyle='--', alpha=0.7)

for i in range(N_modes):
    ax2.plot(s_vals, P_modes[i], color=colors[i], linewidth=2.5, label=f'Mode {i+1}')
ax2.plot(s_vals, P_total, color='black', linestyle='--', linewidth=2, alpha=0.8, label='Total Power')
ax2.set_title(f'Power Dynamics (Exciting Mode 1)', fontsize=14, pad=10)
ax2.set_xlabel(r'Distance $s$, $\mu m$', fontsize=12); ax2.set_ylabel('Power', fontsize=12)
ax2.set_ylim([-0.05, 1.05]); ax2.legend(loc='upper right'); ax2.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout(); plt.savefig("torch_geometry_dynamics.png", dpi=400); plt.show()

# 2. Анализ потерь
fig_loss, (ax_l1, ax_l2) = plt.subplots(2, 1, figsize=(10, 8), dpi=200, sharex=True)
ax_l1.plot(s_vals, loss_db, color='darkred', linewidth=3)
ax_l1.set_ylabel('Cumulative Loss [dB]', fontsize=12)
ax_l1.set_title(f'Insertion Loss Profile (Total: {loss_db[-1]:.3f} dB)', fontsize=14)
ax_l1.grid(True, linestyle='--', alpha=0.7)

ax_l2.fill_between(s_vals, loss_rate, color='coral', alpha=0.4, label='Simulated Loss Rate')
ax_l2.plot(s_vals, loss_rate, color='orangered', linewidth=2)
ax_l2.plot(s_vals, theoretical_rate, color='darkred', linestyle='--', linewidth=2, alpha=0.8, label='Theoretical Bending Loss')
ax_l2.set_ylabel(r'Loss Rate [dB/$\mu m$]', fontsize=12)
ax_l2.set_xlabel(r'Distance $s$, $\mu m$', fontsize=12)
ax_l2.legend(loc='upper right'); ax_l2.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout(); plt.savefig("torch_loss_analysis.png", dpi=400); plt.show()

# 3. 3D Матрицы
def plot_complex_matrix_3d(ax, matrix, title):
    N = matrix.shape[0]
    x, y = np.meshgrid(np.arange(N), np.arange(N))
    x, y = x.flatten(), y.flatten()
    z = np.zeros_like(x)
    magnitudes = np.abs(matrix).flatten()
    phases = np.angle(matrix).flatten()
    cols = cm.hsv((phases + np.pi) / (2 * np.pi))
    dx = dy = 0.6
    ax.bar3d(x - dx/2, y - dy/2, z, dx, dy, magnitudes, color=cols, shade=True, alpha=0.9)
    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel('Input', labelpad=10); ax.set_ylabel('Output', labelpad=10)
    ax.set_zlabel(r'$|U_{ij}|$', labelpad=10)
    ax.set_xticks(np.arange(N)); ax.set_yticks(np.arange(N))
    ax.view_init(elev=30, azim=45)

fig_3d = plt.figure(figsize=(14, 6), dpi=200)
ax1_3d = fig_3d.add_subplot(121, projection='3d')
plot_complex_matrix_3d(ax1_3d, U_target_np, r'Target Unitary $U_{target}$ (DFT)')
ax2_3d = fig_3d.add_subplot(122, projection='3d')
plot_complex_matrix_3d(ax2_3d, M_final, r'Experimental Unitary $M_{sim}$')

sm = cm.ScalarMappable(cmap=cm.hsv, norm=plt.Normalize(vmin=0, vmax=2*np.pi))
cbar_ax = fig_3d.add_axes([0.92, 0.25, 0.015, 0.5])
cbar = fig_3d.colorbar(sm, cax=cbar_ax)
cbar.set_label(r'Phase $\arg(U_{ij})$', fontsize=12)
cbar.set_ticks([0, np.pi/2, np.pi, 3*np.pi/2, 2 * np.pi])
cbar.set_ticklabels([r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
plt.subplots_adjust(wspace=0.1); plt.savefig('torch_unitary_3d.png', dpi=400); plt.show()

# 4. График Кривизны
fig_curv, ax_c = plt.subplots(figsize=(10, 5), dpi=200)
ax_c.plot(s_vals, kappa_final, color='purple', linewidth=2.5, label=r'Optimized $\kappa(s)$')
ax_c.fill_between(s_vals, kappa_final, 0, where=(kappa_final > 0), color='purple', alpha=0.2)
ax_c.fill_between(s_vals, kappa_final, 0, where=(kappa_final <= 0), color='indigo', alpha=0.2)
ax_c.axhline(0, color='black', linewidth=1, linestyle='--')
ax_c.set_title('Optimized Curvature Profile', fontsize=15, pad=10)
ax_c.set_xlabel(r'Distance $s$, $\mu m$', fontsize=12)
ax_c.set_ylabel(r'Curvature $\kappa(s)$, $\mu m^{-1}$', fontsize=12)
ax_c.grid(True, linestyle='--', alpha=0.7); ax_c.legend()
plt.tight_layout(); plt.savefig("torch_curvature.png", dpi=400); plt.show()