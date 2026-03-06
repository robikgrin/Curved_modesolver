import torch
import torch.nn as nn
import numpy as np
from torchdiffeq import odeint

class WaveguideODE(nn.Module):
    """Класс, описывающий производные dU/ds и dPhi/ds для интегратора torchdiffeq"""
    def __init__(self, solver, k_ctrl, L):
        super().__init__()
        self.solver = solver
        self.k_ctrl = k_ctrl
        self.L = L

        self.segments = len(k_ctrl) - 3
        # Матрица базиса кубического B-сплайна
        self.M = torch.tensor([
            [-1,  3, -3,  1],
            [ 3, -6,  3,  0],
            [-3,  0,  3,  0],
            [ 1,  4,  1,  0]
        ], dtype=torch.float64, device=solver.device) / 6.0

    def evaluate_spline(self, s):
        """Вычисляет кривизну и ее производные в конкретной физической точке s (скаляр)"""
        # Переводим физическую длину в глобальный параметр t
        t_global = (s / self.L) * self.segments
        t_global = torch.clamp(t_global, 0.0, self.segments - 1e-5)
        
        idx = torch.floor(t_global).long()
        t = t_global - idx
        
        # Выбираем 4 контрольные точки для текущего локального куска
        P_window = torch.stack([
            self.k_ctrl[idx], 
            self.k_ctrl[idx+1], 
            self.k_ctrl[idx+2], 
            self.k_ctrl[idx+3]
        ])
        
        T = torch.stack([t**3, t**2, t, torch.ones_like(t)])
        dT = torch.stack([3*t**2, 2*t, torch.ones_like(t), torch.zeros_like(t)])
        d2T = torch.stack([6*t, torch.tensor(2.0, dtype=torch.float64, device=self.solver.device), torch.zeros_like(t), torch.zeros_like(t)])
        
        M_P = torch.matmul(self.M, P_window)
        
        k = torch.dot(T, M_P)
        dk_dt = torch.dot(dT, M_P)
        d2k_dt2 = torch.dot(d2T, M_P)
        
        # Пересчет производных (chain rule)
        dt_ds = self.segments / self.L
        dk_ds = dk_dt * dt_ds
        d2k_ds2 = d2k_dt2 * (dt_ds**2)
        
        return k, dk_ds, d2k_ds2

    def forward(self, s, state):
        """Правая часть диффура: возвращает dU/ds и dPhi/ds"""
        U, Phi = state
        
        k, dk, d2k = self.evaluate_spline(s)
        b, db, i1, i2, j0, j1 = self.solver.interpolate_cache(k)
        
        beta_m = b.unsqueeze(1)
        
        # Сборка матриц
        G = 2j * torch.diag(b) + 2 * dk * i1 - dk * j0
        Q = 1j * dk * torch.diag(db) + (dk**2) * i2 + (d2k + 2j * beta_m * dk) * i1 - (dk**2) * j1 - 1j * beta_m * dk * j0
        
        C_raw = torch.linalg.solve(G, -Q)
        
        K = 0.5 * (C_raw - C_raw.mH)
        loss_diag = torch.diag(torch.abs(torch.imag(b)))
        K_lossy = K - loss_diag
        
        # Фазовый синхронизм
        phase_mat = torch.exp(1j * (Phi.unsqueeze(0) - Phi.unsqueeze(1)))
        C_eff = K_lossy * phase_mat
        
        dU_ds = torch.matmul(C_eff, U)
        dPhi_ds = torch.real(b)
        
        return (dU_ds, dPhi_ds)


class TorchModePropagator(nn.Module):
    def __init__(self, cache_file, device='cpu'):
        super().__init__()
        self.device = device
        
        # ЗАГРУЗКА КЭША
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

    def interpolate_cache(self, k):
        # Жесткая защита от выхода за пределы кэша
        k_clamped = torch.clamp(k, self.k_min + 1e-5, self.k_max - 1e-5)
        idx_float = (k_clamped - self.k_min) / self.dk
        idx = torch.floor(idx_float).long()
        w = idx_float - idx
        
        def lerp(tensor):
            return tensor[idx] * (1 - w) + tensor[idx+1] * w

        return lerp(self.beta), lerp(self.db_dk), lerp(self.I1), lerp(self.I2), lerp(self.J0), lerp(self.J1)

    def forward(self, L, k_ctrl):
        ode_func = WaveguideODE(self, k_ctrl, L)
        
        # Начальные условия
        U0 = torch.eye(self.N, dtype=torch.complex128, device=self.device)
        Phi0 = torch.zeros(self.N, dtype=torch.float64, device=self.device)
        
        # Точки, в которых нам нужен результат (начало и конец волновода)
        s_span = torch.tensor([0.0, L], dtype=torch.float64, device=self.device)
        
        # ЗАПУСК АДАПТИВНОГО ИНТЕГРАТОРА DOPRI5 (Адаптивный Рунге-Кутта)
        Y_final = odeint(ode_func, (U0, Phi0), s_span, method='dopri5', rtol=1e-4, atol=1e-5)
        
        # Извлекаем состояние в точке L (последний элемент тензора)
        U_end = Y_final[0][-1]
        Phi_end = Y_final[1][-1]
        
        # Добавляем динамическую фазу
        phase_out = torch.exp(1j * Phi_end)
        U_physical = phase_out.unsqueeze(1) * U_end
        
        if self.training:
            return U_physical
        else:
            # Восстанавливаем плотный вектор кривизны чисто для красивых графиков
            num_points = 400
            segments = len(k_ctrl) - 3
            t_global = torch.linspace(0, segments - 1e-5, num_points, dtype=torch.float64, device=self.device)
            idx = torch.floor(t_global).long()
            t = (t_global - idx).unsqueeze(1)
            
            P_window = torch.stack([k_ctrl[idx], k_ctrl[idx+1], k_ctrl[idx+2], k_ctrl[idx+3]], dim=1)
            T = torch.cat([t**3, t**2, t, torch.ones_like(t)], dim=1)
            M = ode_func.M
            
            M_P = torch.matmul(M, P_window.unsqueeze(-1)).squeeze(-1)
            k_dense = torch.sum(T * M_P, dim=1)
            
            return U_physical, k_dense