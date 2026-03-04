import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

class TorchModePropagator(nn.Module):
    def __init__(self, cache_file, num_steps=1000, device='cpu'):
        super().__init__()
        self.device = device
        self.num_steps = num_steps
        
        # 1. ЗАГРУЗКА КЭША
        data = np.load(cache_file)
        self.k_array = torch.tensor(data['kappa_array'], dtype=torch.float64, device=device)
        self.dk = self.k_array[1] - self.k_array[0]
        self.k_min = self.k_array[0]
        self.k_max = self.k_array[-1]
        self.N = data['beta'].shape[1] 
        
        self.beta = torch.tensor(data['beta'], dtype=torch.complex128, device=device)
        self.db_dk = torch.tensor(data['db_dk'], dtype=torch.complex128, device=device)
        self.I1 = torch.tensor(data['I1'], dtype=torch.complex128, device=device)
        self.I2 = torch.tensor(data['I2'], dtype=torch.complex128, device=device)
        self.J0 = torch.tensor(data['J0'], dtype=torch.complex128, device=device)
        self.J1 = torch.tensor(data['J1'], dtype=torch.complex128, device=device)

    def bezier_curve(self, P, num_points):
        """
        Генерирует гладкую кривую Безье и ее производные аналитически.
        P - тензор контрольных точек [N_ctrl]
        """
        n = len(P) - 1
        t = torch.linspace(0, 1, num_points, dtype=torch.float64, device=self.device)
        
        # Вычисляем полиномы Бернштейна
        B = torch.zeros((num_points, n + 1), dtype=torch.float64, device=self.device)
        dB = torch.zeros((num_points, n), dtype=torch.float64, device=self.device)
        d2B = torch.zeros((num_points, n - 1), dtype=torch.float64, device=self.device)
        
        import math
        for i in range(n + 1):
            coef = math.comb(n, i)
            B[:, i] = coef * (t ** i) * ((1 - t) ** (n - i))
            
        for i in range(n):
            coef = n * math.comb(n - 1, i)
            dB[:, i] = coef * (t ** i) * ((1 - t) ** (n - 1 - i))
            
        for i in range(n - 1):
            coef = n * (n - 1) * math.comb(n - 2, i)
            d2B[:, i] = coef * (t ** i) * ((1 - t) ** (n - 2 - i))
            
        k = torch.matmul(B, P)
        dk_dt = torch.matmul(dB, P[1:] - P[:-1])
        d2k_dt2 = torch.matmul(d2B, P[2:] - 2*P[1:-1] + P[:-2])
        
        return k, dk_dt, d2k_dt2

    def interpolate_cache(self, k):
        # Мягкий Clamp, который не убивает градиенты мгновенно
        k_clamped = k - F.relu(k - self.k_max + 1e-5) + F.relu(self.k_min + 1e-5 - k)
        
        idx = (k_clamped - self.k_min) / self.dk
        idx_floor = torch.floor(idx).long()
        idx_floor = torch.clamp(idx_floor, 0, len(self.k_array) - 2)
        
        w = (idx - idx_floor)
        
        def lerp(tensor):
            if tensor.dim() == 2:
                return tensor[idx_floor] * (1 - w.unsqueeze(-1)) + tensor[idx_floor+1] * w.unsqueeze(-1)
            else:
                return tensor[idx_floor] * (1 - w.unsqueeze(-1).unsqueeze(-1)) + tensor[idx_floor+1] * w.unsqueeze(-1).unsqueeze(-1)

        return lerp(self.beta), lerp(self.db_dk), lerp(self.I1), lerp(self.I2), lerp(self.J0), lerp(self.J1)

    def get_C_matrix(self, k, dk, d2k, Phi):
        b, db, i1, i2, j0, j1 = self.interpolate_cache(k)
        beta_m = b.unsqueeze(1) # shape: [Batch, N, 1]
        
        # Растягиваем скаляры dk, d2k до нужной размерности
        dk = dk.unsqueeze(-1).unsqueeze(-1)
        d2k = d2k.unsqueeze(-1).unsqueeze(-1)
        
        G = 2j * torch.diag_embed(b) + 2 * dk * i1 - dk * j0
        Q = 1j * dk * torch.diag_embed(db) + (dk**2) * i2 + (d2k + 2j * beta_m * dk) * i1 - (dk**2) * j1 - 1j * beta_m * dk * j0
        
        # Пакетное (Batched) решение системы: G * C_raw = -Q
        C_raw = torch.linalg.solve(G, -Q)
        
        # Симметризация
        K = 0.5 * (C_raw - C_raw.mH)
        loss_diag = torch.diag_embed(torch.abs(torch.imag(b)))
        K_lossy = K - loss_diag
        
        phase_mat = torch.exp(1j * (Phi.unsqueeze(-1) - Phi.unsqueeze(-2)))
        return K_lossy * phase_mat, torch.real(b)

    def forward(self, L, P_ctrl):
        ds = L / self.num_steps
        
        # Генерируем гладкие кривые Безье
        k_vals, dk_dt, d2k_dt2 = self.bezier_curve(P_ctrl, self.num_steps)
        
        # Пересчитываем производные по t в производные по s (ds = L * dt)
        dk_vals = dk_dt / L
        d2k_vals = d2k_dt2 / (L**2)
        
        U = torch.eye(self.N, dtype=torch.complex128, device=self.device)
        Phi = torch.zeros(self.N, dtype=torch.float64, device=self.device)
        
        for i in range(self.num_steps):
            k = k_vals[i].unsqueeze(0)
            dk = dk_vals[i].unsqueeze(0)
            d2k = d2k_vals[i].unsqueeze(0)
            
            # RK4
            C1, dP1 = self.get_C_matrix(k, dk, d2k, Phi.unsqueeze(0))
            C1, dP1 = C1[0], dP1[0]
            k1_U = torch.matmul(C1, U)
            
            C2, dP2 = self.get_C_matrix(k, dk, d2k, (Phi + dP1 * ds/2).unsqueeze(0))
            C2, dP2 = C2[0], dP2[0]
            k2_U = torch.matmul(C2, U + k1_U * ds/2)
            
            C3, dP3 = self.get_C_matrix(k, dk, d2k, (Phi + dP2 * ds/2).unsqueeze(0))
            C3, dP3 = C3[0], dP3[0]
            k3_U = torch.matmul(C3, U + k2_U * ds/2)
            
            C4, dP4 = self.get_C_matrix(k, dk, d2k, (Phi + dP3 * ds).unsqueeze(0))
            C4, dP4 = C4[0], dP4[0]
            k4_U = torch.matmul(C4, U + k3_U * ds)
            
            U = U + (ds / 6.0) * (k1_U + 2*k2_U + 2*k3_U + k4_U)
            Phi = Phi + (ds / 6.0) * (dP1 + 2*dP2 + 2*dP3 + dP4)
            
        return U