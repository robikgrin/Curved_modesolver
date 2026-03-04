
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from scipy.sparse import csc_matrix, diags, bmat
from scipy.sparse.linalg import eigs


class rect_WG:
    r'''Class of rectangular profile waveguide

    Attributes
    ----------
    ``wavelength`` : float
            Wavelength of input electromagnetic wave
    
    ``n_clad`` : float
            Refractive index of cladding
    
    ``n_core`` : float
            Refractive index of core
    
    ``d_xi`` : float
            Simulation step in ``xi`` direction
    
    ``d_eta`` : float
            Simulation step in ``eta`` drection

    ``W`` : float
            Width of Si core
    
    ``H`` : float
            Height of Si core

    ``delta_l`` : float
            Distance to the left border of simulation
    
    ``delta_r`` : float
            Distance to the right border of simulation
    
    ``delta_u`` : float
            Distance to the upper border of simulation
    
    ``delta_d`` : float
            Distance to the down border of simulation
    
    ``kappa`` : float
            Curvature value (``kappa = 0`` as default)

    Methods
    ----------
    ``draw_structure()`` :
            Draw refractive index profile of defined waveguide structure
    
    ``draw_permittivity_profile()`` :
            Draw permittivity profiles ``ex, ey, ez``, calculated by index averaging technique
    
    ``get_grid_info()`` : 
            Get information about number of grids in xi direction ``N``, eta direction ``Q``, upper-left grid coordinate of Si core ``n_left, q_up``, lower-rigth grid of Si core ``n_rigth, q_down``
    
    ``get_PML(NPML)`` : 
            Setting PML

    ``draw_PML(proj)`` : 
            Visualization of PML components for defined projection: imaginary and real part 
    
    ``FDE(num)`` :
            Finite-difference eigensolver, which calculates eigenmodes and eigenvalues of curved waveguide
    
    ``get_field(mode_num)`` :
            Calculating the field projections for each projection for defined mode
    
    ``get_E_field(mode_num)`` : 
            Electric density prjections normalizator for defined mode
    
    ``draw_field(mode_num, scale)`` : 
            Draw a density field projections for both electric and magnetic density profiles for orthonormal curvilinear axes ``xi``, ``eta`` and ``s`` in defined scale
    
    ``draw_E_filed(mode_num, scale)`` : 
            Draw a electric density field projections for orthonormal curvilinear axes ``xi``, ``eta`` and ``s`` in defined scale
    
    ``get_overlap(num, lap_type)`` :
            Overlap calculations between all modes for defined overlap formula
    '''
    def __init__(self, wavelength=1.55E-6, n_clad=1.444, n_core=3.4755, d_xi=0.02E-6, d_eta=0.02E-6, W=2E-6, H = 0.22E-6, delta_l=2E-6, delta_r=2E-6, delta_u=2E-6, delta_d=2E-6, kappa=0):
        r'''Sets the initial parameters of waveguide structure for eigenmode calculation

        Parameters
        -----------
        ``wavelength`` : float
                Wavelength of input electromagnetic wave. ``wavelength = 1.55E-6`` as default
        
        ``n_clad`` : float
                Refractive index of cladding ( ``n_clad = 1.444`` as default (SiO2 cladding) )

        ``n_core`` : float
                Refractive index of core ( ``n_core = 3.4755`` as default (Si core) )

        ``d_xi`` : float
                Simulation step in ``xi`` direction ( ``d_xi = 0.02E-6`` as default)

        ``d_xi`` : float
                Simulation step in ``eta`` drection ( ``d_eta = 0.02E-6`` as default)

        ``W`` : float
                Width of Si core ( ``W = 2E-6`` as default)

        ``H`` : float
                Height of Si core (``H = 0.22E-6`` as default)

        ``delta_l`` : float
                Distance to the left border of simulation (``delta_l = 2E-6`` as default)

        ``delta_r`` : float
                Distance to the right border of simulation (``delta_r = 2E-6`` as default)

        ``delta_u`` : float
                Distance to the upper border of simulation (``delta_u = 2E-6`` as default)

        ``delta_d`` : float
                Distance to the down border of simulation (``delta_d = 2E-6`` as default)

        ``kappa`` : float
                Curvature value (``kappa = 0`` as default)
        '''
        self.k0 = 2*np.pi/wavelength
        self.n_clad = n_clad
        self.n_core = n_core
        self.d_xi = d_xi
        self.d_eta = d_eta
        self.d_X = d_xi*self.k0
        self.d_Y = d_eta*self.k0
        self.W = W
        self.H = H
        self.kappa = kappa
        
        if self.kappa >= 0:
            self.delta_r = delta_r
            self.delta_l = delta_l
        else:
            self.delta_r = delta_l
            self.delta_l = delta_r
        
        self.delta_u = delta_u
        self.delta_d = delta_d

        x_size = delta_r + delta_l + W  #size of full location in \xi direction
        y_size = H + delta_u + delta_d #size of full location in \eta direction

        N = int(x_size/d_xi) - 1  #number of grids in \xi direction (horizontal direction)
        Q = int(y_size/d_eta) - 1 #number of grids in \eta direction (vertical direction)

        ### Silicon structure ###'

        n_l = int(self.delta_l/d_xi) #left-up grid (\xi number)
        q_u = int(self.delta_u/d_eta) #lef-up grid (\eta number)
        n_r = int((W+self.delta_l)/d_xi) #right-down grid (\xi number)
        q_d = int((self.delta_u + H)/d_eta) #right-down grid (\eta number)

        self.n_l = n_l
        self.n_r = n_r
        self.q_u = q_u
        self.q_d = q_d

        ### the n matrix ###
        n_index = np.full((Q, N), n_clad, dtype=float)
        n_index[q_u:q_d, n_l:n_r] = n_core
        
        self.n_index = n_index
        self.N = N
        self.Q = Q

        new_ind = np.full((Q+2, N+2), n_clad**2)
        new_ind[1:Q+1, 1:N+1] = n_index**2

        e_rx = np.copy(new_ind)
        for q in range(1, Q+1):
            for n in range(1, N+1):
                e_rx[q, n] = 0.5*(new_ind[q, n] + new_ind[q-1, n])

        e_ry = np.copy(new_ind)
        for q in range(1, Q+1):
            for n in range(1, N+1):
                e_ry[q, n] = 0.5*(new_ind[q, n] + new_ind[q, n-1])

        e_rz = np.copy(new_ind)
        for q in range(1, Q+1):
            for n in range(1, N+1):
                e_rz[q, n] = 0.25*(new_ind[q, n] + new_ind[q-1, n-1]+ new_ind[q-1, n] + new_ind[q, n-1])

        self.ex = np.copy(e_rx[1:Q+1, 1:N+1])
        self.ey = np.copy(e_ry[1:Q+1, 1:N+1])
        self.ez = np.copy(e_rz[1:Q+1, 1:N+1])

        _xi_ = np.linspace(-W/2 - delta_l, W/2 + delta_r, N)
        self.t_xi = 1 + kappa * _xi_
    
    def set_width(self, W:float):
        r'''Set width of core value

        Parameters
        ----------
        ``W`` : float
                New width value
        '''
        self.W = W

    def draw_structure(self):
        r'''Draw refractive index profile of defined waveguide structure'''
        x_ticks = np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N)
        y_ticks = np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q)
        X, Y = np.meshgrid(x_ticks * 1e6, y_ticks * 1e6)

        fig, ax = plt.subplots()
        c = ax.pcolormesh(X, Y, self.n_index, cmap='crest_r', shading='auto')
        fig.colorbar(c, ax=ax)
        
        ax.set_title(r'The area of reaserch (in each cell value of $n^{q}_n$)')
        ax.set_xlabel(r'$\xi$ coordinate, $\mu m$')
        ax.set_ylabel(r'$\eta$ coordinate, $\mu m$')
        ax.set_aspect('equal')
        plt.show()

    def draw_permittivity_profile(self):
        r'''Draw permittivity profiles ``ex, ey, ez``, calculated by index averaging technique'''
        x_ticks = np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N)
        y_ticks = np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q)
        X, Y = np.meshgrid(x_ticks * 1e6, y_ticks * 1e6)

        titles = [r'$\epsilon_{\xi}$ profile in each cell', 
                  r'$\epsilon_{\eta}$ profile in each cell', 
                  r'$\epsilon_{s}$ profile in each cell']
        data = [self.ex, self.ey, self.ez]

        for d, title in zip(data, titles):
            fig, ax = plt.subplots()
            c = ax.pcolormesh(X, Y, d, cmap='crest_r', shading='auto')
            fig.colorbar(c, ax=ax)
            ax.set_title(title)
            ax.set_xlabel(r'$\xi$ coordinate, $\mu m$')
            ax.set_ylabel(r'$\eta$ coordinate, $\mu m$')
            ax.set_aspect('equal')
            plt.show()
    
    def get_grid_info(self):
        r'''Get information about number of grids in xi direction ``N``, eta direction ``Q``, upper-left grid coordinate of Si core ``n_left, q_up``, lower-rigth grid of Si core ``n_rigth, q_down``
        '''
        n_l = self.n_l
        q_u = self.q_u
        n_r = self.n_core
        q_d = self.q_d

        print(f'Number of grids in xi direction: N = {self.N} \nNumber of grids in eta direction: Q = {self.Q}')

        print(f"The (n,q) value for left upper grid of Si: ({n_l}, {q_u})")
        print(f"The (n,q) value for the right lower grid of Si: ({n_r}, {q_d})")

    def get_PML(self, NPML):
        r'''Setting PML 

        Parameters
        ------------
        ``NPML`` : list or ndarray
                The size of PML layer in each direcion ``[x_left, x_right, y_down, y_up]``
        '''
        # x axis
        sx = np.ones((self.Q, self.N), dtype=complex)
        
        # Левый PML
        n_left = np.arange(NPML[0])
        profile_left = (1 + 3 * ((NPML[0] - n_left) / NPML[0])**3) * \
                       (1 + 1j * self.n_clad * (np.sin(np.pi * (NPML[0] - n_left) / (2 * NPML[0])))**2)
        # Для строк NumPy сам умеет "растягивать" вектор по вертикали
        sx[:, :NPML[0]] = profile_left
        
        # Правый PML
        n_right = np.arange(self.N - NPML[1], self.N)
        profile_right = (1 + 3 * ((n_right - (self.N - NPML[1])) / NPML[1])**3) * \
                        (1 + 1j * self.n_clad * (np.sin(np.pi * (n_right - (NPML[1] + self.N)) / (2 * NPML[1])))**2)
        sx[:, -NPML[1]:] = profile_right

        # y axis
        sy = np.ones((self.Q, self.N), dtype=complex)
        
        # Верхний PML (или нижний, в зависимости от твоей системы координат)
        q_up = np.arange(self.Q - NPML[2], self.Q)
        profile_up = (1 + 3 * ((q_up - (self.Q - NPML[2])) / NPML[2])**3) * \
                     (1 + 1j * self.n_clad * (np.sin(np.pi * (q_up - (self.Q - NPML[2])) / (2 * NPML[2])))**2)
        # Превращаем 1D вектор в вектор-столбец с помощью [:, None], чтобы размножить по горизонтали
        sy[-NPML[2]:, :] = profile_up[:, None]

        # Нижний PML
        q_down = np.arange(NPML[3])
        profile_down = (1 + 3 * ((NPML[3] - q_down) / NPML[3])**3) * \
                       (1 + 1j * self.n_clad * (np.sin(np.pi * (NPML[3] - q_down) / (2 * NPML[3])))**2)
        sy[:NPML[3], :] = profile_down[:, None]

        self.sx = sx
        self.sy = sy
    
    def draw_PML(self, proj:str):
        r'''Visualization of PML components for defined projection: imaginary and real part'''
        x_ticks = np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N)
        y_ticks = np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q)
        X, Y = np.meshgrid(x_ticks * 1e6, y_ticks * 1e6)

        if proj == 'x':
            data = self.sx
            title_real = r'Real part of $s_{x}(x)$'
            title_imag = r'Imaginary part of $s_{x}(x)$'
        elif proj == 'y':
            data = self.sy
            title_real = r'Real part of $s_{y}(y)$'
            title_imag = r'Imaginary part of $s_{y}(y)$'
        else:
            print('ERROR: Wrong projection (possible projections are x or y)')
            return

        for d, title in zip([np.real(data), np.imag(data)], [title_real, title_imag]):
            fig, ax = plt.subplots()
            c = ax.pcolormesh(X, Y, d, cmap='crest_r', shading='auto')
            fig.colorbar(c, ax=ax)
            ax.set_title(title)
            ax.set_xlabel(r'$\xi$ coordinate, $\mu m$')
            ax.set_ylabel(r'$\eta$ coordinate, $\mu m$')
            ax.set_aspect('equal')
            plt.show()
    
    @staticmethod
    def U_xx(_N: int, _Q: int, D_X: float):
        size = _N * _Q
        
        # Создаем массивы для диагоналей
        main_diag = np.full(size, -1.0 / D_X)
        up_diag = np.full(size - 1, 1.0 / D_X)
        
        # Обнуляем элементы на границах строк, чтобы производная 
        # по X не "перепрыгивала" на следующий ряд по Y
        up_diag[_N - 1 :: _N] = 0.0 
        
        # Твое краевое условие для последнего элемента
        main_diag[-1] = 1.0 / D_X 

        # Собираем разреженную матрицу одним вызовом
        return diags([main_diag, up_diag], offsets=[0, 1], format='csc', dtype=complex)

    @staticmethod
    def U_yy(_N: int, _Q: int, D_Y: float):
        size = _N * _Q
        
        # Главная диагональ
        main_diag = np.full(size, -1.0 / D_Y)
        # Верхняя диагональ (сдвиг на _N, так как берем соседнюю ячейку по Y)
        up_diag = np.full(size - _N, 1.0 / D_Y)
        
        return diags([main_diag, up_diag], offsets=[0, _N], format='csc', dtype=complex)

    @staticmethod
    def V_xx(_N: int, _Q: int, D_X: float):
        size = _N * _Q
        
        # Главная диагональ
        main_diag = np.full(size, 1.0 / D_X)
        # Нижняя диагональ (сдвиг -1)
        low_diag = np.full(size - 1, -1.0 / D_X)
        
        # Обнуляем элементы на границах строк, чтобы производная не "перепрыгивала" по Y
        low_diag[_N - 1 :: _N] = 0.0 
        
        return diags([low_diag, main_diag], offsets=[-1, 0], format='csc', dtype=complex)

    @staticmethod
    def V_yy(_N: int, _Q: int, D_Y: float):
        size = _N * _Q
        
        # Главная диагональ
        main_diag = np.full(size, 1.0 / D_Y)
        # Нижняя диагональ (сдвиг на -_N)
        low_diag = np.full(size - _N, -1.0 / D_Y)
        
        return diags([low_diag, main_diag], offsets=[-_N, 0], format='csc', dtype=complex)

    @staticmethod
    def S_xi(_N: int, _Q: int, _sx, _sy, _t_xi):
        # Автоматическое умножение 2D массива на 1D вектор и вытягивание в линию
        main_diag = ((_sy / _sx) * _t_xi).flatten()
        return diags([main_diag], offsets=[0], format='csc', dtype=complex)

    @staticmethod
    def S_eta(_N: int, _Q: int, _sx, _sy, _t_xi):
        main_diag = ((_sx / _sy) * _t_xi).flatten()
        return diags([main_diag], offsets=[0], format='csc', dtype=complex)

    @staticmethod
    def inv_S_s(_N: int, _Q: int, _sx, _sy):
        main_diag = (1.0 / (_sx * _sy)).flatten()
        return diags([main_diag], offsets=[0], format='csc', dtype=complex)
    
    @staticmethod
    def e_xi(_N: int, _Q: int, e_x):
        return diags([e_x.flatten()], offsets=[0], format='csc', dtype=complex)

    @staticmethod
    def e_eta(_N: int, _Q: int, e_y):
        return diags([e_y.flatten()], offsets=[0], format='csc', dtype=complex)

    @staticmethod
    def e_s_inverse(_N: int, _Q: int, e_z):
        return diags([(1.0 / e_z).flatten()], offsets=[0], format='csc', dtype=complex)

    def T_xi(self):
        main_diag = (self.ex * (self.sy / self.sx) * self.t_xi).flatten()
        return diags([main_diag], offsets=[0], format='csc', dtype=complex)

    def T_eta(self):
        main_diag = (self.ey * (self.sx / self.sy) * self.t_xi).flatten()
        return diags([main_diag], offsets=[0], format='csc', dtype=complex)

    def inv_T_s(self):
        main_diag = (1.0 / (self.ez * self.sx * self.sy)).flatten()
        return diags([main_diag], offsets=[0], format='csc', dtype=complex)

    def A_xx(self):
        # Главная диагональ: дублируем массив -t_xi / d_X ровно Q раз
        main_diag = np.tile(-self.t_xi / self.d_X, self.Q)
        main_diag[-1] = self.t_xi[-1] / self.d_X  # Твое краевое условие
        
        # Верхняя диагональ: берем тот же массив, но со сдвигом +1 (срез [1:])
        up_diag = np.tile(self.t_xi / self.d_X, self.Q)[1:]
        # Обнуляем переходы между строками
        up_diag[self.N - 1 :: self.N] = 0.0
        
        return diags([main_diag, up_diag], offsets=[0, 1], format='csc', dtype=complex)

    def A_yy(self):
        # Главная диагональ
        main_diag = np.tile(-self.t_xi / self.d_Y, self.Q)
        
        # Верхняя диагональ (сдвиг +N, поэтому дублируем на 1 раз меньше: Q-1)
        up_diag = np.tile(self.t_xi / self.d_Y, self.Q - 1)
        
        return diags([main_diag, up_diag], offsets=[0, self.N], format='csc', dtype=complex)

    def C_xx(self):
        # Главная диагональ
        main_diag = np.tile(self.t_xi / self.d_X, self.Q)
        
        # Нижняя диагональ: сдвиг -1, поэтому срез [:-1]
        low_diag = np.tile(-self.t_xi / self.d_X, self.Q)[:-1]
        # Обнуляем переходы между строками
        low_diag[self.N - 1 :: self.N] = 0.0
        
        return diags([low_diag, main_diag], offsets=[-1, 0], format='csc', dtype=complex)

    def C_yy(self):
        # Главная диагональ
        main_diag = np.tile(self.t_xi / self.d_Y, self.Q)
        
        # Нижняя диагональ (сдвиг -N, дублируем на 1 раз меньше: Q-1)
        low_diag = np.tile(-self.t_xi / self.d_Y, self.Q - 1)
        
        return diags([low_diag, main_diag], offsets=[-self.N, 0], format='csc', dtype=complex)

    def FDE(self, num: int, NPML: list = None, v0 = None):
        r'''Finite-difference eigensolver, which calculates eigenmodes and eigenvalues of curved waveguide
        
        Parameters
        --------
        num : int 
                number of calculating modes
        NPML : list
                The size of PML layer in each direcion [x_left, x_right, y_down, y_up]
        '''
        # Если NPML не передали при вызове, ставим твои старые значения по умолчанию
        if NPML is None:
            NPML = [150, 150, 60, 60]

        U_x, U_y, V_x, V_y = self.U_xx(self.N, self.Q, self.d_X), self.U_yy(self.N, self.Q, self.d_Y), self.V_xx(self.N, self.Q, self.d_X), self.V_yy(self.N, self.Q, self.d_Y)

        A_x = self.A_xx()
        A_y = self.A_yy()
        C_x = self.C_xx()
        C_y = self.C_yy()
        
        self.get_PML(NPML)
        
        S_xi, S_eta, inv_S_s = self.S_xi(self.N, self.Q, self.sx, self.sy, self.t_xi), self.S_eta(self.N, self.Q, self.sx, self.sy, self.t_xi), self.inv_S_s(self.N, self.Q, self.sx, self.sy)
        T_xi, T_eta, inv_T_s = self.T_xi(), self.T_eta(), self.inv_T_s()
            
        P_xx = -A_x @ inv_T_s @ V_y @ C_x @ inv_S_s @ U_y + (S_eta + A_x @ inv_T_s @ V_x) @ (T_xi + C_y @ inv_S_s @ U_y)
        P_xy = A_x @ inv_T_s @ V_y @ (T_eta + C_x @ inv_S_s @ U_x) - (S_eta + A_x @ inv_T_s @ V_x) @ C_y @ inv_S_s @ U_x
        P_yx = -(S_xi + A_y @ inv_T_s @ V_y) @ C_x @ inv_S_s @ U_y + A_y @ inv_T_s @ V_x @ (T_xi + C_y @ inv_S_s @ U_y)
        P_yy = (T_eta + C_x @ inv_S_s @ U_x) @ (S_xi + A_y @ inv_T_s @ V_y) - A_y @ inv_T_s @ V_x @ C_y @ inv_S_s @ U_x
        
        P = bmat([[P_xx, P_xy], 
                  [P_yx, P_yy]], format='csc')

        ### EIGENVALUES ###
        target_neff_sq = self.n_core**2
        if v0 is None:
            vals, vecs = eigs(P, k=num, sigma=target_neff_sq, which='LM')
        else:
            vals, vecs = eigs(P, k=num, sigma=target_neff_sq, which='LM', v0=v0)
        n_eff = np.sqrt(vals)
        self.n_eff = n_eff
        self.vecs = vecs

    def get_field(self, mode_num:int):
        r'''Calculating the field projections for each projection for defined mode

        Parameters
        -----------
        mode : int
                Mode number

        Return:
        -----------------
        ``E_x, E_y, E_z, H_x, H_y, H_z`` - density profile output in matrix ``Q * N`` form

        ``h_x, h_y`` - magnetic density projections for ``xi`` and ``eta`` axes
        '''
        U_x, U_y, V_x, V_y = self.U_xx(self.N, self.Q, self.d_X), self.U_yy(self.N, self.Q, self.d_Y), self.V_xx(self.N, self.Q, self.d_X), self.V_yy(self.N, self.Q, self.d_Y)
        
        C_x, C_y = self.C_xx(), self.C_yy()
        T_xi, T_eta, inv_T_s = self.T_xi(), self.T_eta(),  self.inv_T_s()
        inv_S_s = self.inv_S_s(self.N, self.Q, self.sx, self.sy)

        h_x = (1/self.n_eff[mode_num - 1]) * (C_x @ inv_S_s @ U_y @ self.vecs[:self.N*self.Q, mode_num-1] - (T_eta + C_x @ inv_S_s @ U_x) @ self.vecs[self.N*self.Q : 2*self.N*self.Q, mode_num-1])
        h_y = (1/self.n_eff[mode_num - 1]) * ((T_xi + C_y @ inv_S_s @ U_y) @ self.vecs[:self.N*self.Q, mode_num-1] - C_y @ inv_S_s @ U_x @ self.vecs[self.N*self.Q : 2*self.N*self.Q, mode_num-1])
        h_z = 1j * inv_S_s @ (-U_y @ self.vecs[:self.N*self.Q, mode_num-1] + U_x @ self.vecs[self.N*self.Q : 2*self.N*self.Q, mode_num-1])

        e_z = 1j * inv_T_s @ (V_y @ h_x - V_x @ h_y)
        #in matrix form
        E_x = np.zeros((self.Q, self.N), dtype=complex)
        E_y = np.zeros((self.Q, self.N), dtype=complex)
        E_z = np.zeros((self.Q, self.N), dtype=complex)

        H_x = np.zeros((self.Q, self.N), dtype=complex)
        H_y = np.zeros((self.Q, self.N), dtype=complex)
        H_z = np.zeros((self.Q, self.N), dtype=complex)

        for i in range(self.Q):
            E_x[i, :] = self.vecs[i*self.N:(i+1)*self.N, mode_num - 1]
            E_y[i, :] = self.vecs[self.N*self.Q + i*self.N : self.N*self.Q + (i+1)*self.N, mode_num-1]
            E_z[i, :] = e_z[i*self.N:(i+1)*self.N]

            H_x[i, :] = h_x[i*self.N:(i+1)*self.N]
            H_y[i, :] = h_y[i*self.N:(i+1)*self.N]
            H_z[i, :] = h_z[i*self.N:(i+1)*self.N]
            
        return E_x, E_y, E_z, H_x, H_y, H_z, h_x, h_y

    def get_E_field(self, mode_num:int):
        r'''Electric density prjections normalizator

        Parameters
        ----------
        ``mode`` : int 
                mode number 
        ``scale`` : str 
                scaling format for data visualization (`log` or `norm`)

        Return:
        -----------------
        ``E_field`` - density profile output with three components in one vector
        '''
        size = self.N * self.Q

        # Генерируем операторы (В будущем их лучше вынести в __init__, чтобы не пересчитывать каждый раз)
        U_x, U_y = self.U_xx(self.N, self.Q, self.d_X), self.U_yy(self.N, self.Q, self.d_Y)
        V_x, V_y = self.V_xx(self.N, self.Q, self.d_X), self.V_yy(self.N, self.Q, self.d_Y)
        C_x, C_y = self.C_xx(), self.C_yy()
        T_xi, T_eta, inv_T_s = self.T_xi(), self.T_eta(),  self.inv_T_s()
        inv_S_s = self.inv_S_s(self.N, self.Q, self.sx, self.sy)

        # Вытаскиваем сразу ВСЕ нужные моды (форма матриц: size x mode_num)
        vecs_x = self.vecs[:size, :mode_num]
        vecs_y = self.vecs[size : 2*size, :mode_num]

        # Инвертированный n_eff для всей пачки мод сразу
        n_inv = 1.0 / self.n_eff[:mode_num]

        # Считаем h_x и h_y матрично для всех мод сразу. 
        # Скобки расставлены так, чтобы вычислять "справа налево" — это в разы быстрее!
        h_x = (C_x @ (inv_S_s @ (U_y @ vecs_x)) - (T_eta + C_x @ inv_S_s @ U_x) @ vecs_y) * n_inv
        h_y = ((T_xi + C_y @ inv_S_s @ U_y) @ vecs_x - C_y @ (inv_S_s @ (U_x @ vecs_y))) * n_inv

        # e_z для всех мод
        e_z = 1j * inv_T_s @ (V_y @ h_x - V_x @ h_y)

        # Считаем норму вдоль столбцов (axis=0) сразу для всех мод
        norm_sq = np.sum(np.abs(vecs_x)**2, axis=0) + \
                  np.sum(np.abs(vecs_y)**2, axis=0) + \
                  np.sum(np.abs(e_z)**2, axis=0)
        norm = np.sqrt(norm_sq * self.d_X * self.d_Y)
        
        # Нормируем компоненты
        vec_norm_x = vecs_x / norm
        vec_norm_y = vecs_y / norm
        vec_norm_z = e_z / norm

        # Собираем в итоговую матрицу (форма: 3*size x mode_num)
        E_field = np.vstack((vec_norm_x, vec_norm_y, vec_norm_z))
        
        return E_field

    def draw_field(self, mode:int, scale:str):
        r'''Draw a density field projections for both electric and magnetic density profiles...'''
        from matplotlib.colors import LogNorm, Normalize
        
        x_ticks = np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N)
        y_ticks = np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q)
        X, Y = np.meshgrid(x_ticks * 1e6, y_ticks * 1e6)

        E_x, E_y, E_z, H_x, H_y, H_z = self.get_field(mode_num=mode)[:6]
        
        fields = [[E_x, E_y, E_z], 
                  [H_x, H_y, H_z]]
        titles = [[r'$E_{\xi}$', r'$E_{\eta}$', r'$E_{s}$'],
                  [r'$H_{\xi}$', r'$H_{\eta}$', r'$H_{s}$']]

        fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(24, 10),
                               gridspec_kw={'width_ratios': [1, 1, 1], 'height_ratios': [5, 5],
                                            'wspace': 0.2, 'hspace': 0.4})

        # Твоя оригинальная шапка
        fig.suptitle(r'Field projections for ' + str(mode) +  r' mode : $n_{eff}$ =' + str(np.round(self.n_eff[mode-1], 3)) + r", $\kappa = $" + str(self.kappa) + r" $\mu m^{-1}$", fontsize=20)

        for i in range(2):
            for j in range(3):
                ax = axs[i, j]
                
                # Создаем новый нормализатор для КАЖДОГО графика внутри цикла!
                norm = LogNorm() if scale == 'log' else Normalize()
                
                # Отрисовка поля
                c = ax.pcolormesh(X, Y, np.abs(fields[i][j]), cmap='jet', shading='auto', norm=norm)
                
                # Твои подписи
                ax.set_title(titles[i][j], fontsize=15)
                ax.set_xlabel(r'$\xi$ coordinate, $\mu m$')
                ax.set_ylabel(r'$\eta$ coordinate, $\mu m$')
                ax.set_aspect('equal')
                fig.colorbar(c, ax=ax, fraction=0.046, pad=0.04)

        plt.show()

    def draw_E_field(self, mode:int, scale:str):
        r'''Draw a electric density field projections for orthonormal curvilinear axes...'''
        from matplotlib.colors import LogNorm, Normalize
        
        x_ticks = np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N)
        y_ticks = np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q)
        X, Y = np.meshgrid(x_ticks * 1e6, y_ticks * 1e6)

        E = self.get_E_field(mode_num=mode)
        E_x = E[:self.N*self.Q, mode-1].reshape(self.Q, self.N)
        E_y = E[self.N*self.Q : 2*self.N*self.Q, mode-1].reshape(self.Q, self.N)
        E_z = E[2*self.N*self.Q : 3*self.N*self.Q, mode-1].reshape(self.Q, self.N)

        fields = [E_x, E_y, E_z]
        titles = [r'$E_{\xi}$', r'$E_{\eta}$', r'$E_{s}$']

        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(24, 10),
                               gridspec_kw={'width_ratios': [1, 1, 1], 'height_ratios': [5],
                                            'wspace': 0.2, 'hspace': 0.4})

        # Твоя оригинальная шапка
        fig.suptitle(r'Field projections for ' + str(mode) +  r' mode : $n_{eff}$ =' + str(np.round(self.n_eff[mode-1], 3)) + r", $\kappa = $" + str(self.kappa) + r" $\mu m^{-1}$", fontsize=20)

        for j in range(3):
            ax = axs[j]
            data_to_plot = np.abs(fields[j])**2 if scale == 'norm' else np.abs(fields[j])
            
            # Создаем новый нормализатор для КАЖДОГО графика внутри цикла!
            norm = LogNorm() if scale == 'log' else Normalize()
            
            c = ax.pcolormesh(X, Y, data_to_plot, cmap='jet', shading='auto', norm=norm)
            ax.set_title(titles[j], fontsize=15)
            ax.set_xlabel(r'$\xi$ coordinate, $\mu m$')
            ax.set_ylabel(r'$\eta$ coordinate, $\mu m$')
            ax.set_aspect('equal')
            fig.colorbar(c, ax=ax, fraction=0.046, pad=0.04)

        plt.show()
    
    @staticmethod
    def lum_summa(E_xi, E_eta, H_xi, H_eta):
        return np.dot(E_xi, np.conj(H_eta)) - np.dot(np.conj(H_xi), E_eta)

    @staticmethod
    def my_summa(x, y):
        return np.dot(x, np.conj(y))
    
    def get_overlap(self, num:int, lap_type: str):
        r'''Overlap calculations between all modes for defined overlap formula

        Parameters
        ----------
        ``num`` : int
                Number of calculated modes from ``FDE(num)`` function
        
        ``lap_type`` : str
                Overlap formula defining: ``lumerical`` - ANSYS Lumerical formula, ``classic`` - classical overlap formula

        Return
        ------
        ``res`` : pd.Dataframe
                Overlap dataframe for all modes
        '''
        G_mat = np.zeros((num, num), dtype=float)
        if lap_type == 'lumerical':
            for m in range(1, num+1):
                for n in range(1, m+1):
                    H_xi_m, H_eta_m = self.get_field(m)[6:8]
                    H_xi_n, H_eta_n = self.get_field(n)[6:8]
                    E_xi_m, E_eta_m = self.vecs[:self.N*self.Q, m-1], self.vecs[self.N*self.Q : 2*self.N*self.Q, m-1]
                    E_xi_n, E_eta_n = self.vecs[:self.N*self.Q, n-1], self.vecs[self.N*self.Q : 2*self.N*self.Q, n-1]

                    upper = np.real(self.lum_summa(E_xi_m, E_eta_m, H_xi_n, H_eta_n) * self.lum_summa(E_xi_n, E_eta_n, H_xi_m, H_eta_m)/self.lum_summa(E_xi_m, E_eta_m, H_xi_m, H_eta_m))
                    down = np.real(self.lum_summa(E_xi_n, E_eta_n, H_xi_n, H_eta_n))

                    G_mat[m-1, n-1] = np.abs(upper/down)/self.k0
                    G_mat[n-1, m-1] = G_mat[m-1, n-1]
        
        elif lap_type == 'classic':
            for m in range(1, num+1):
                for n in range(1, m+1):
                    E_xi_m, E_eta_m = self.vecs[:self.N*self.Q, m-1], self.vecs[self.N*self.Q : 2*self.N*self.Q, m-1]
                    E_xi_n, E_eta_n = self.vecs[:self.N*self.Q, n-1], self.vecs[self.N*self.Q : 2*self.N*self.Q, n-1]

                    upper = np.real(self.my_summa(E_xi_m, E_xi_n) + self.my_summa(E_eta_m, E_eta_n))
                    down = np.sqrt(self.my_summa(E_xi_m, E_xi_m) + self.my_summa(E_eta_m, E_eta_m)) * np.sqrt(self.my_summa(E_xi_n, E_xi_n) + self.my_summa(E_eta_n, E_eta_n))

                    G_mat[m-1, n-1] = np.abs(upper/down)/self.k0
                    G_mat[n-1, m-1] = G_mat[m-1, n-1]    
        #plotting
        mat = np.copy(G_mat)
        
        for i in range(num):
            mat[i,i] = 1
        
        res = pd.DataFrame(mat, index=[i for i in range(1, num+1)], columns=[i for i in range(1, num+1)])

        ax = plt.axes()
        sns.heatmap(res, ax = ax, linewidths=.7, linecolor='white', cmap='Spectral_r', cbar=True, norm=LogNorm(), robust = True, fmt='',
                    annot_kws={'color': 'black'})
        ax.set_title(r'Overlap for all calculated eigenmodes profiles (lumerical)')
        ax.set_xlabel(r'Mode number')
        ax.set_ylabel(r'Mode number')
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')

        plt.show()

        return res