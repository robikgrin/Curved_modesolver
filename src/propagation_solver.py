import numpy as np
import curved_modesolver as cms
from scipy.integrate import solve_ivp
from tqdm import tqdm

class ModePropagator:
    r'''
    Solver for mode propagation in a curved multimode waveguide.
    Calculates overlap integrals and coupling coefficients.
    '''
    def __init__(self, wg_params: dict, num_modes: int, NPML: list = None, d_kappa: float = 50):
        r'''
        Parameters
        ----------
        `wg_params` : dict
            Словарь с параметрами для инициализации rect_WG (без kappa)
        `num_modes` : int
            Количество мод для расчета
        `NPML` : list
            Параметры PML слоя
        `d_kappa` : float
            Шаг для численного дифференцирования по кривизне (d\kappa)
        '''
        self.params = wg_params
        self.num_modes = num_modes
        self.NPML = NPML if NPML is not None else [100, 100, 100, 100]
        self.dk = d_kappa
        self.v0 = None
    def _solve_for_kappa(self, kappa):
        r"""Вспомогательный метод: решает FDE для конкретного \kappa с горячим стартом"""
        wg = cms.rect_WG(**self.params, kappa=kappa)
        
        wg.FDE(self.num_modes, self.NPML, v0=self.v0)

        if hasattr(wg, 'evecs'):
            self.v0_guess = np.real(wg.evecs[:, 0])
        
        beta = wg.k0 * wg.n_eff[:self.num_modes]
        U = wg.get_E_field(self.num_modes)
        
        return wg, beta, U

    def get_derivatives(self, kappa):
        r"""Численное дифференцирование по кривизне методом центральных разностей"""
        # 1. Решаем для центральной точки (kappa)
        wg_0, b_0, U_0 = self._solve_for_kappa(kappa)
        
        # 2. Решаем для (kappa - d_kappa)
        _, b_m, U_m = self._solve_for_kappa(kappa - self.dk)
        
        # 3. Решаем для (kappa + d_kappa)
        _, b_p, U_p = self._solve_for_kappa(kappa + self.dk)

        # Первая производная (центральная разность O(h^2))
        db_dk = (b_p - b_m) / (2 * self.dk)
        dU_dk = (U_p - U_m) / (2 * self.dk)
        
        # Вторая производная
        d2U_dk2 = (U_p - 2*U_0 + U_m) / (self.dk**2)

        return wg_0, b_0, db_dk, U_0, dU_dk, d2U_dk2

    def compute_integrals(self, kappa):
        r'''
        Вычисляет матрицы интегралов перекрытия I^(1), I^(2), J^(0), J^(1).
        input: `kappa` - значение кривизны для расчета
        return: Возвращает словарик с матрицами размера (`num_modes` x `num_modes`).
        '''
        wg, beta, db_dk, U, dU_dk, d2U_dk2 = self.get_derivatives(kappa)
        
        N, Q = wg.N, wg.Q
        dS = wg.d_X * wg.d_Y

        # Строим плоский массив координаты \xi, вытянутый для N*Q ячеек
        x_ticks = np.linspace(-wg.W/2 - wg.delta_l, wg.W/2 + wg.delta_r, N)
        xi_2d = np.tile(x_ticks, (Q, 1))
        xi_flat = xi_2d.flatten()
        

        xi_full = np.tile(xi_flat, 3) 
        
        geom_factor = xi_full / (1 + kappa * xi_full)

        I1 = np.zeros((self.num_modes, self.num_modes), dtype=complex)
        I2 = np.zeros((self.num_modes, self.num_modes), dtype=complex)
        J0 = np.zeros((self.num_modes, self.num_modes), dtype=complex)
        J1 = np.zeros((self.num_modes, self.num_modes), dtype=complex)

        for n in range(self.num_modes):
            U_n_conj = np.conj(U[:, n])
            
            for m in range(self.num_modes):
                I1[n, m] = np.sum(U_n_conj * dU_dk[:, m]) * dS
                
                I2[n, m] = np.sum(U_n_conj * d2U_dk2[:, m]) * dS
                
                J0[n, m] = np.sum(geom_factor * U_n_conj * U[:, m]) * dS
                
                J1[n, m] = np.sum(geom_factor * U_n_conj * dU_dk[:, m]) * dS

        return {
            'beta': beta,
            'db_dk': db_dk,
            'I1': I1,
            'I2': I2,
            'J0': J0,
            'J1': J1
        }
    
    def auto_detect_modes(self, max_search=5):
        r'''Автоматически определяет количество физически направляемых мод.'''
        k0 = 2 * np.pi / self.params['wavelength']
        n_clad = self.params['n_clad']
        
        # Временно устанавливаем большое количество мод, чтобы найти их все
        self.num_modes = max_search
        
        # Считаем сырые матрицы для прямого волновода (kappa = 0)
        # Используем твой внутренний метод (наверняка он называется compute_integrals или аналогично)
        res = self.compute_integrals(0.0) 
        betas = np.real(res['beta'])
        neffs = betas / k0
        
        # Считаем, сколько мод удовлетворяют условию полного внутреннего отражения
        physical_count = np.sum(neffs > n_clad)
        
        print("\n" + "="*40)
        print(" ФИЗИЧЕСКИЙ АНАЛИЗ ВОЛНОВОДА (kappa = 0)")
        print("="*40)
        for i, neff in enumerate(neffs):
            if neff > n_clad:
                print(f" Мода {i+1}: n_eff = {neff:.4f}  [+] Направляемая")
            else:
                print(f" Мода {i+1}: n_eff = {neff:.4f}  [-] Вытекающая в PML")
                
        print(f"\nУстановлено рабочее количество мод: {physical_count}")
        print("="*40 + "\n")
        
        self.num_modes = physical_count
        return physical_count
    
    def calculate_and_save_cache(self, kappa_array, filename="wg_cache.npz"):
        r'''Считает интегралы для сетки значений kappa и сохраняет их в бинарный файл.'''
        print(f"Начинаем тяжелый расчет матриц для {len(kappa_array)} значений kappa...")
        data = {key: [] for key in ['beta', 'db_dk', 'I1', 'I2', 'J0', 'J1']}
        
        for k in tqdm(kappa_array, desc="Расчет мод и интегралов", unit="точек"):
            res = self.compute_integrals(k)
            for key in data.keys():
                data[key].append(res[key])
                
        save_dict = {'kappa_array': kappa_array}
        for key in data.keys():
            save_dict[key] = np.array(data[key])
            
        np.savez_compressed(filename, **save_dict)
        print(f"\nРасчет завершен. Матрицы сохранены в: {filename}")

    def load_cache_and_build_splines(self, filename="wg_cache.npz"):
        r'''Загружает предрасчитанные матрицы и строит ОДИН быстрый мега-сплайн.'''
        print(f"Загрузка данных из {filename}...")
        data = np.load(filename)
        kappa_array = data['kappa_array']
        
        N = self.num_modes
        packed_data = []
        
        # Сплющиваем все матрицы в один длинный вектор для каждого значения kappa
        for i in range(len(kappa_array)):
            vec = np.concatenate([
                data['beta'][i].flatten(),
                data['db_dk'][i].flatten(),
                data['I1'][i].flatten(),
                data['I2'][i].flatten(),
                data['J0'][i].flatten(),
                data['J1'][i].flatten()
            ])
            packed_data.append(vec)
            
        packed_data = np.array(packed_data)
        
        # Строим ОДИН кубический сплайн для всего вектора сразу
        from scipy.interpolate import CubicSpline
        self.mega_spline = CubicSpline(kappa_array, packed_data, axis=0)
        print("Быстрый Мега-сплайн построен!")

    def solve_dynamics(self, s_span, kappa_func, dkappa_func, d2kappa_func, c0, fast_mode=False):
        r'''Решает систему ОДУ, принудительно сохраняя энергию (Hermitian Coupled Mode Theory).'''
        N = self.num_modes

        # Нам нужны только амплитуды и фазы
        Phi0 = np.zeros(N, dtype=float)
        y0 = np.concatenate([c0, Phi0]).astype(complex)

        def dY_ds(s, Y):
            c = Y[0:N]
            Phi = np.real(Y[N:2*N]) # Фаза строго вещественная

            k = kappa_func(s)
            dk = dkappa_func(s)
            d2k = d2kappa_func(s)

            vec = self.mega_spline(k)
            beta = vec[0:N]
            db_dk = vec[N:2*N]
            I1 = vec[2*N : 2*N + N**2].reshape((N, N))
            I2 = vec[2*N + N**2 : 2*N + 2*N**2].reshape((N, N))
            J0 = vec[2*N + 2*N**2 : 2*N + 3*N**2].reshape((N, N))
            J1 = vec[2*N + 3*N**2 : 2*N + 4*N**2].reshape((N, N))

            beta_m = beta[np.newaxis, :]
            
            # Строим сырые матрицы G и Q
            G = 2j * np.diag(beta) + 2 * dk * I1 - dk * J0
            Q = 1j * dk * np.diag(db_dk) + (dk**2) * I2 + (d2k + 2j * beta_m * dk) * I1 - (dk**2) * J1 - 1j * beta_m * dk * J0

            # 1. Получаем "сырую" матрицу связи
            C_raw = np.linalg.solve(G, -Q)

            # 2. ПРИНУДИТЕЛЬНЫЙ ЗАКОН СОХРАНЕНИЯ ЭНЕРГИИ (Симметризация)
            # Мы отсекаем весь численный шум, оставляя только анти-эрмитовую часть (перекачку)
            K = 0.5 * (C_raw - np.conj(C_raw.T))

            # 3. Добавляем строгое физическое затухание (от PML)
            # Вычитаем модуль мнимой части бета, чтобы гарантировать потерю, а не усиление
            loss_diag = np.diag(np.abs(np.imag(beta)))
            K_lossy = K - loss_diag

            # 4. Учитываем фазовый синхронизм (Phase matching)
            phase_mat = np.exp(1j * (Phi[np.newaxis, :] - Phi[:, np.newaxis]))
            C_eff = K_lossy * phase_mat

            # Вычисляем производные
            dc = C_eff @ c
            dPhi = np.real(beta)

            return np.concatenate([dc, dPhi])

        if fast_mode:
            sol = solve_ivp(dY_ds, s_span, y0, method='RK45', dense_output=False, rtol=1e-3, atol=1e-6)
        else:
            sol = solve_ivp(dY_ds, s_span, y0, method='RK45', dense_output=True, rtol=1e-6, atol=1e-8)
        
        return sol