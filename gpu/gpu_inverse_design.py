import torch
import torch.nn as nn
import torch.nn.functional as F
import matplotlib.pyplot as plt
import numpy as np
from gpu_solver import TorchModePropagator
import time

# --- НАСТРОЙКИ ---
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Вычисления запущены на: {device}")

# Имя твоего файла кэша (скопированного в эту папку!)
CACHE_FILE = "wg_cache_micro.npz" 
L = 70.0
num_steps = 300 # Плотная сетка интегратора

# Инициализируем наш PyTorch солвер
solver = TorchModePropagator(CACHE_FILE, num_steps=num_steps, device=device)

# --- ПАРАМЕТРЫ ВОЛНОВОДА КАК "ВЕСА НЕЙРОСЕТИ" ---
# Даем 15 контрольных точек кривизны. Требуем вычисления градиентов (requires_grad=True)
N_ctrl = 15
k_ctrl = nn.Parameter(torch.zeros(N_ctrl, dtype=torch.float64, device=device))

# Оптимизатор Adam (классика глубокого обучения)
optimizer = torch.optim.Adam([k_ctrl], lr=0.01)

# --- ЦЕЛЕВАЯ МАТРИЦА ---
# Мы хотим делитель Адамара (50/50 по мощности с правильной фазой)
# Для 2-х мод матрица выглядит так:
Target_M = (1 / np.sqrt(2)) * torch.tensor([
    [1.0,  1j],
    [1j,  1.0]
], dtype=torch.complex128, device=device)

print("\nЗапуск Градиентного Обратного Проектирования (Adjoint Method)...")
start_time = time.time()

# --- ЦИКЛ ОБУЧЕНИЯ (TRAINING LOOP) ---
epochs = 200
loss_history = []

for epoch in range(epochs):
    optimizer.zero_grad()
    
    # 1. Гладкая интерполяция (Upsampling) контрольных точек в плотный массив для ОДУ
    # PyTorch сам дифференцирует этот шаг! Добавляем нули на концы (clamped)
    k_full = torch.cat([torch.tensor([0.0], device=device), k_ctrl, torch.tensor([0.0], device=device)])
    k_dense = F.interpolate(k_full.view(1, 1, -1), size=num_steps, mode='linear', align_corners=True).squeeze()
    
    # 2. Прямой проход (Forward pass) через солвер ОДУ
    # На выходе сразу вся матрица передачи!
    M_sim = solver(L, k_dense)
    
    # 3. Функция потерь (Loss)
   # Размерность нашей матрицы (N_modes)
    N = Target_M.shape[0]
    
    # Вычисляем след произведения эрмитово-сопряженной цели на нашу матрицу: Tr(U^dagger * M)
    trace_val = torch.trace(torch.matmul(Target_M.mH, M_sim))
    
    # Считаем фиделити (от 0.0 до 1.0)
    fidelity = (torch.abs(trace_val)**2) / (N**2)
    
    # Ошибка - это отклонение фиделити от идеальной единицы
    matrix_error = 1.0 - fidelity
    
    # Регуляризация (Total Variation), чтобы волновод был гладким
    smoothness_penalty = 0.05 * torch.sum((k_ctrl[1:] - k_ctrl[:-1])**2)
    
    loss = matrix_error + smoothness_penalty
    
    if (epoch + 1) % 20 == 0:
        print(f"Итерация {epoch+1}/{epochs} | Loss: {loss.item():.4f} | Fidelity: {fidelity.item():.4f}")
    
print(f"Оптимизация завершена за {time.time() - start_time:.2f} секунд!")

# --- ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ ---
k_final = k_dense.detach().cpu().numpy()
s_vals = np.linspace(0, L, num_steps)

# Получаем итоговую матрицу на CPU
M_final = M_sim.detach().cpu().numpy()

print("\nИтоговая Матрица Передачи (Модули в квадрате - мощности):")
print(np.round(np.abs(M_final)**2, 3))

plt.figure(figsize=(12, 4))
plt.subplot(1, 2, 1)
plt.plot(loss_history, 'm-', linewidth=2)
plt.title('Кривая падения ошибки (Loss)')
plt.xlabel('Эпоха оптимизации (Epoch)')
plt.ylabel('Loss')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(s_vals, k_final, 'b-', linewidth=2)
plt.title('Оптимальная функция кривизны $\kappa(s)$')
plt.xlabel('Координата $s$, мкм')
plt.ylabel('Кривизна $\kappa$, мкм$^{-1}$')
plt.grid(True)
plt.tight_layout()
plt.show()