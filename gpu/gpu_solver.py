import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

class TorchModePropagator(nn.Module):
    def __init__(self, cache_file, num_steps=200, device='cpu'):
        super().__init__()
        self.device = device
        self.num_steps = num_steps
        
        # 1. ЗАГРУЗКА КЭША В ТЕНЗОРЫ
        data = np.load(cache_file)
        self.k_array = torch.tensor(data['kappa_array'], dtype=torch.float64, device=device)
        self.dk = self.k_array[1] - self.k_array[0]
        self.k_min = self.k_array[0]
        self.k_max = self.k_array[-1]
        
        self.N = data['beta'].shape[1] # Количество мод
        
        # Переносим все матрицы на GPU/CPU в комплексном формате
        self.beta = torch.tensor(data['beta'], dtype=torch.complex128, device=device)
        self.db_dk = torch.tensor(data['db_dk'], dtype=torch.complex128, device=device)
        self.I1 = torch.tensor(data['I1'], dtype=torch.complex128, device=device)
        self.I2 = torch.tensor(data['I2'], dtype=torch.complex128, device=device)
        self.J0 = torch.tensor(data['J0'], dtype=torch.complex128, device=device)
        self.J1 = torch.tensor(data['J1'], dtype=torch.complex128, device=device)

    def interpolate_cache(self, k):
        """Быстрая дифференцируемая линейная интерполяция кэша на GPU"""
        # Защита от выхода за границы кэша
        k_clamped = torch.clamp(k, self.k_min + 1e-6, self.k_max - 1e-6)
        
        # Вычисляем дробные индексы для интерполяции
        idx = (k_clamped - self.k_min) / self.dk
        idx_floor = torch.floor(idx).long()
        weight = (idx - idx_floor).unsqueeze(-1).unsqueeze(-1) if self.I1.dim() > 2 else (idx - idx_floor).unsqueeze(-1)
        
        def lerp_matrix(tensor, i, w):
            if tensor.dim() == 2: # Для векторов (beta, db_dk)
                return tensor[i] * (1 - w.squeeze(-1)) + tensor[i+1] * w.squeeze(-1)
            else: # Для матриц NxN
                return tensor[i] * (1 - w) + tensor[i+1] * w

        b = lerp_matrix(self.beta, idx_floor, weight)
        db = lerp_matrix(self.db_dk, idx_floor, weight)
        i1 = lerp_matrix(self.I1, idx_floor, weight)
        i2 = lerp_matrix(self.I2, idx_floor, weight)
        j0 = lerp_matrix(self.J0, idx_floor, weight)
        j1 = lerp_matrix(self.J1, idx_floor, weight)
        
        return b, db, i1, i2, j0, j1

    def get_C_matrix(self, k, dk_ds, d2k_ds2, Phi):
        """Сборка консервативной матрицы связи (с поддержкой Autograd)"""
        b, db, i1, i2, j0, j1 = self.interpolate_cache(k)
        
        beta_m = b.unsqueeze(0)
        
        G = 2j * torch.diag(b) + 2 * dk_ds * i1 - dk_ds * j0
        Q = 1j * dk_ds * torch.diag(db) + (dk_ds**2) * i2 + (d2k_ds2 + 2j * beta_m * dk_ds) * i1 - (dk_ds**2) * j1 - 1j * beta_m * dk_ds * j0
        
        # Решаем систему G * C_raw = -Q
        C_raw = torch.linalg.solve(G, -Q)
        
        # Принудительное сохранение энергии (анти-эрмитова часть)
        K = 0.5 * (C_raw - C_raw.mH)
        
        # Добавляем физические потери
        loss_diag = torch.diag(torch.abs(torch.imag(b)))
        K_lossy = K - loss_diag
        
        # Учет фазы
        phase_mat = torch.exp(1j * (Phi.unsqueeze(1) - Phi.unsqueeze(0)))
        C_eff = K_lossy * phase_mat
        
        return C_eff, torch.real(b)

    def forward(self, L, kappa_dense):
        """Прогон диффура RK4. Выдает полную матрицу передачи."""
        ds = L / self.num_steps
        
        # Считаем производные кривизны (конечные разности на GPU)
        dk_ds = torch.gradient(kappa_dense, spacing=(ds,))[0]
        d2k_ds2 = torch.gradient(dk_ds, spacing=(ds,))[0]
        
        # Стартуем с ЕДИНИЧНОЙ матрицы. Каждый столбец - это отдельный порт входа!
        U = torch.eye(self.N, dtype=torch.complex128, device=self.device)
        Phi = torch.zeros(self.N, dtype=torch.float64, device=self.device)
        
        # Метод Рунге-Кутты 4-го порядка в чистом тензорном виде
        for i in range(self.num_steps):
            k_val = kappa_dense[i]
            dk_val = dk_ds[i]
            d2k_val = d2k_ds2[i]
            
            # Шаг 1
            C1, dP1 = self.get_C_matrix(k_val, dk_val, d2k_val, Phi)
            k1_U = torch.matmul(C1, U)
            
            # Шаг 2
            C2, dP2 = self.get_C_matrix(k_val, dk_val, d2k_val, Phi + dP1 * ds/2)
            k2_U = torch.matmul(C2, U + k1_U * ds/2)
            
            # Шаг 3
            C3, dP3 = self.get_C_matrix(k_val, dk_val, d2k_val, Phi + dP2 * ds/2)
            k3_U = torch.matmul(C3, U + k2_U * ds/2)
            
            # Шаг 4
            C4, dP4 = self.get_C_matrix(k_val, dk_val, d2k_val, Phi + dP3 * ds)
            k4_U = torch.matmul(C4, U + k3_U * ds)
            
            # Обновление состояния
            U = U + (ds / 6.0) * (k1_U + 2*k2_U + 2*k3_U + k4_U)
            Phi = Phi + (ds / 6.0) * (dP1 + 2*dP2 + 2*dP3 + dP4)
            
        return U