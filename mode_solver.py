
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from scipy.sparse import csc_matrix
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
    
    ``d_xi`` : float
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

        for q in range(q_u, q_d):
            for n in range(n_l, n_r):
                n_index[q, n] = n_core
        
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
        r'''Draw refractive index profile of defined waveguide structure
        '''
        x_ticks = np.round(10**6 * np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N), 1)
        y_ticks = np.round(10**6 * np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q), 1)

        df = pd.DataFrame(data=self.n_index, index=y_ticks, columns=x_ticks)
        ax = plt.axes()

        sns.heatmap(df, ax=ax, xticklabels=self.N//4, yticklabels=self.Q//7, cmap='crest_r')
        ax.set_title(r'The area of reaserch (in each cell value of $n^{q}_n$)')
        ax.set_xlabel(r'$\xi$ coordinate, $\mu m$')
        ax.set_ylabel(r'$\eta$ coordinate, $\mu m$')
        plt.yticks(rotation=0) 
        plt.show()
  
    def draw_permittivity_profile(self):
        r'''Draw permittivity profiles ``ex, ey, ez``, calculated by index averaging technique
        '''
        x_ticks = np.round(10**6 * np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N), 1)
        y_ticks = np.round(10**6 * np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q), 1)

        df_x= pd.DataFrame(data=self.ex, index=y_ticks, columns=x_ticks)
        df_y= pd.DataFrame(data=self.ey, index=y_ticks, columns=x_ticks)
        df_z= pd.DataFrame(data=self.ez, index=y_ticks, columns=x_ticks)
        
        ax_1 = plt.axes()
        sns.heatmap(df_x, ax=ax_1, xticklabels=self.N//4, yticklabels=self.Q//7, cmap='crest_r')
        ax_1.set_title(r'$\epsilon_{\xi}$ profile in each cell')
        ax_1.set_xlabel(r'$\xi$ coordinate, $\mu m$')
        ax_1.set_ylabel(r'$\eta$ coordinate, $\mu m$')
        plt.yticks(rotation=0) 
        plt.show()

        ax_2 = plt.axes()
        sns.heatmap(df_y, ax=ax_2, xticklabels=self.N//4, yticklabels=self.Q//7, cmap='crest_r')
        ax_2.set_title(r'$\epsilon_{\eta}$ profile in each cell')
        ax_2.set_xlabel(r'$\xi$ coordinate, $\mu m$')
        ax_2.set_ylabel(r'$\eta$ coordinate, $\mu m$')
        plt.yticks(rotation=0) 
        plt.show()

        ax_3 = plt.axes()
        sns.heatmap(df_z, ax=ax_3, xticklabels=self.N//4, yticklabels=self.Q//7, cmap='crest_r')
        ax_3.set_title(r'$\epsilon_{s}$ profile in each cell')
        ax_3.set_xlabel(r'$\xi$ coordinate, $\mu m$')
        ax_3.set_ylabel(r'$\eta$ coordinate, $\mu m$')
        plt.yticks(rotation=0) 
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
        #x axis
        sx = np.ones((self.Q, self.N), dtype=complex)
        for n in range(NPML[0]):
            sx[:,n] = (1 + 3*((NPML[0] - n)/NPML[0])**3) * (1 + 1j*self.n_clad*(np.sin(np.pi*(NPML[0] - n)/(2*NPML[0])))**2)
        
        for n in range(self.N - NPML[1], self.N):
            sx[:,n] = (1 + 3*((n - (self.N - NPML[1]))/NPML[1])**3) * (1 + 1j*self.n_clad*(np.sin(np.pi*(n - (NPML[1] + self.N))/(2*NPML[1])))**2)

        #y axis
        sy = np.ones((self.Q, self.N), dtype=complex)
        for q in range(self.Q - NPML[2], self.Q):
            sy[q, :] = (1 + 3*((q - (self.Q - NPML[2]))/NPML[2])**3) * (1 + 1j*self.n_clad*(np.sin(np.pi*(q - (self.Q - NPML[2]))/(2*NPML[2])))**2)

        for q in range(NPML[3]):
            sy[q, :] = (1 + 3*((NPML[3] - q)/NPML[3])**3) * (1 + 1j*self.n_clad*(np.sin(np.pi*(NPML[3] - q)/(2*NPML[3])))**2)

        self.sx = sx
        self.sy = sy
    
    def draw_PML(self, proj:str):
        r'''Visualization of PML components for defined projection: imaginary and real part 
        
        Parameters
        --------
        ``proj`` : str
                projection (``x`` or ``y``)
        '''
        x_ticks = np.round(10**6 * np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N), 1)
        y_ticks = np.round(10**6 * np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q), 1)

        if proj == 'x':
            df_1 = pd.DataFrame(data=self.sx, index=y_ticks, columns=x_ticks)
            
            ax_1 = plt.axes()
            sns.heatmap(np.real(df_1), ax=ax_1, xticklabels=self.N//4, yticklabels=self.Q//7, cmap='crest_r')
            ax_1.set_title(r'Imaginary part of $s_{x}(x)$')
            ax_1.set_xlabel(r'$\xi$ coordinate, $\mu m$')
            ax_1.set_ylabel(r'$\eta$ coordinate, $\mu m$')
            plt.yticks(rotation=0) 
            plt.show()

            ax_2 = plt.axes()
            sns.heatmap(np.imag(df_1), ax=ax_2, xticklabels=self.N//4, yticklabels=self.Q//7, cmap='crest_r')
            ax_2.set_title(r'Imaginary part of $s_{x}(x)$')
            ax_2.set_xlabel(r'$\xi$ coordinate, $\mu m$')
            ax_2.set_ylabel(r'$\eta$ coordinate, $\mu m$')
            plt.yticks(rotation=0) 
            plt.show()

        elif proj == 'y':
            df_1 = pd.DataFrame(data=self.sy, index=y_ticks, columns=x_ticks)
            
            ax_1 = plt.axes()
            sns.heatmap(np.real(df_1), ax=ax_1, xticklabels=self.N//4, yticklabels=self.Q//7, cmap='crest_r')
            ax_1.set_title(r'Imaginary part of $s_{y}(y)$')
            ax_1.set_xlabel(r'$\xi$ coordinate, $\mu m$')
            ax_1.set_ylabel(r'$\eta$ coordinate, $\mu m$')
            plt.yticks(rotation=0) 
            plt.show()

            ax_2 = plt.axes()
            sns.heatmap(np.imag(df_1), ax=ax_2, xticklabels=self.N//4, yticklabels=self.Q//7, cmap='crest_r')
            ax_2.set_title(r'Imaginary part of $s_{y}(y)$')
            ax_2.set_xlabel(r'$\xi$ coordinate, $\mu m$')
            ax_2.set_ylabel(r'$\eta$ coordinate, $\mu m$')
            plt.yticks(rotation=0) 
            plt.show()
        else:
            print('ERROR: Wrong projection (possible projections are x or y)')
    @staticmethod
    def U_xx(_N: int, _Q: int, D_X: float):
        row_main = np.array([i for i in range(_N*_Q)])
        column_main = np.array([i for i in range(_N*_Q)])
        value_main = np.array([-1/(D_X) for i in range(_N*_Q)])

        row_up = np.array([])
        column_up = np.array([])
        value_up =  np.array([])

        for q in range(_Q):
            temp_row_up = np.array([i for i in range(q*_N, (q+1)*_N-1)])
            temp_column_up = np.array([i for i in range(q*_N + 1, (q+1)*_N)])
            temp_value_up = np.array([1/D_X for n in range(1, _N)])

            value_up = np.append(value_up, temp_value_up)

            row_up = np.append(row_up, temp_row_up)
            column_up = np.append(column_up, temp_column_up)

        value_main[-1] = 1/D_X
        row = np.append(row_main, row_up)
        column = np.append(column_main, column_up)
        value = np.append(value_main, value_up)

        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)

    @staticmethod
    def U_yy(_N: int, _Q: int, D_Y: float):
        row_main = np.array([i for i in range(_N*_Q)])
        column_main = np.array([i for i in range(_N*_Q)])
        value_main = np.array([-1/(D_Y) for i in range(_N*_Q)])

        row_up = np.array([i for i in range(_N*(_Q-1))])
        column_up = np.array([i for i in range(_N, _N*_Q)])
        value_up = np.array([1/(D_Y) for i in range(_N*(_Q-1))])

        row = np.append(row_main, row_up)
        column = np.append(column_main, column_up)
        value = np.append(value_main, value_up)

        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)
    
    @staticmethod
    def V_xx(_N: int, _Q: int, D_X: float):
        row_main = np.array([i for i in range(_N*_Q)])
        column_main = np.array([i for i in range(_N*_Q)])
        value_main = np.array([1/(D_X) for i in range(_N*_Q)])

        row_low = np.array([])
        column_low = np.array([])
        value_low = np.array([])

        for q in range(_Q):
            temp_row_low = np.array([i for i in range(q*_N + 1, (q+1)*_N)])
            temp_column_low = np.array([i for i in range(q*_N, (q+1)*_N-1)])
            temp_value_low = np.array([-1/D_X for n in range(0, _N-1)])

            value_low = np.append(value_low, temp_value_low)

            row_low = np.append(row_low, temp_row_low)
            column_low = np.append(column_low, temp_column_low)

        row = np.append(row_main, row_low)
        column = np.append(column_main, column_low)
        value = np.append(value_main, value_low)

        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)

    @staticmethod
    def V_yy(_N: int, _Q: int, D_Y: float):
        row_main = np.array([i for i in range(_N*_Q)])
        column_main = np.array([i for i in range(_N*_Q)])
        value_main = np.array([1/(D_Y) for i in range(_N*_Q)])

        row_low = np.array([i for i in range(_N, _N*_Q)])
        column_low = np.array([i for i in range(_N*(_Q-1))])
        value_low = np.array([-1/(D_Y) for i in range(_N*(_Q-1))])

        row = np.append(row_main, row_low)
        column = np.append(column_main, column_low)
        value = np.append(value_main, value_low)

        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)

    @staticmethod
    def S_xi(_N: int, _Q: int, _sx, _sy, _t_xi):
        temp = np.zeros((_Q, _N), dtype=complex)
        for q in range(_Q):
            for n in range(_N):
                temp[q,n] = (_sy[q,n]/_sx[q,n]) * _t_xi[n]
        
        row = np.array([i for i in range(_N*_Q)])
        column = np.array([i for i in range(_N*_Q)])
        
        value = np.array([0 for i in range(_N*_Q)], dtype=complex)
        for i in range(_Q):
            value[i*_N: (i+1)*_N] = temp[i, :]
        
        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)
    
    @staticmethod
    def S_eta(_N: int, _Q: int, _sx, _sy, _t_xi):
        temp = np.zeros((_Q, _N), dtype=complex)
        for q in range(_Q):
            for n in range(_N):
                temp[q,n] = (_sx[q,n]/_sy[q,n]) * _t_xi[n]
        
        row = np.array([i for i in range(_N*_Q)])
        column = np.array([i for i in range(_N*_Q)])
        
        value = np.array([0 for i in range(_N*_Q)], dtype=complex)
        for i in range(_Q):
            value[i*_N: (i+1)*_N] = temp[i, :]
        
        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)
    
    @staticmethod
    def inv_S_s(_N: int, _Q: int, _sx, _sy):
        temp = np.zeros((_Q, _N), dtype=complex)
        for q in range(_Q):
            for n in range(_N):
                temp[q,n] = _sx[q,n]*_sy[q,n]
        
        row = np.array([i for i in range(_N*_Q)])
        column = np.array([i for i in range(_N*_Q)])
        
        value = np.array([0 for i in range(_N*_Q)], dtype=complex)
        for i in range(_Q):
            value[i*_N: (i+1)*_N] = 1/temp[i, :]
        
        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)
    
    @staticmethod
    def e_xi(_N: int, _Q: int, e_x):
        row = np.array([i for i in range(_N*_Q)])
        column = np.array([i for i in range(_N*_Q)])
        
        value = np.array([0 for i in range(_N*_Q)], dtype= float)
        for i in range(_Q):
            value[i*_N: (i+1)*_N] = e_x[i, :]

        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)

    @staticmethod
    def e_eta(_N: int, _Q: int, e_y):
        row = np.array([i for i in range(_N*_Q)])
        column = np.array([i for i in range(_N*_Q)])
        
        value = np.array([0 for i in range(_N*_Q)], dtype= float)
        for i in range(_Q):
            value[i*_N : (1+i)*_N] = e_y[i, :]

        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)

    @staticmethod
    def e_s_inverse(_N: int, _Q: int, e_z):
        temp = 1/np.copy(e_z)

        #sns.heatmap(e_rz)
        row = np.array([i for i in range(_N*_Q)])
        column = np.array([i for i in range(_N*_Q)])
        
        value = np.array([0 for i in range(_N*_Q)], dtype= float)
        for i in range(_Q):
            value[i*_N : (1+i)*_N] = temp[i, :]

        return csc_matrix((value, (row, column)), shape=(_N*_Q, _N*_Q), dtype=complex)

    def T_xi(self):
        temp = np.zeros((self.Q, self.N), dtype=complex)
        for q in range(self.Q):
            for n in range(self.N):
                temp[q,n] = self.ex[q,n] * (self.sy[q,n]/self.sx[q,n]) * self.t_xi[n]
        T_temp = np.copy(temp)

        row = np.array([i for i in range(self.N*self.Q)])
        column = np.array([i for i in range(self.N*self.Q)])
        
        value = np.array([0 for i in range(self.N*self.Q)], dtype=complex)
        for i in range(self.Q):
            value[i*self.N: (i+1)*self.N] = T_temp[i, :]
        
        return csc_matrix((value, (row, column)), shape=(self.N*self.Q, self.N*self.Q), dtype=complex)

    def T_eta(self):
        temp = np.zeros((self.Q, self.N), dtype=complex)
        for q in range(self.Q):
            for n in range(self.N):
                temp[q,n] = self.ey[q,n]*(self.sx[q,n]/self.sy[q,n]) * self.t_xi[n]
        
        T_temp = np.copy(temp)

        row = np.array([i for i in range(self.N*self.Q)])
        column = np.array([i for i in range(self.N*self.Q)])
        
        value = np.array([0 for i in range(self.N*self.Q)], dtype=complex)
        for i in range(self.Q):
            value[i*self.N: (i+1)*self.N] = T_temp[i, :]
        
        return csc_matrix((value, (row, column)), shape=(self.N*self.Q, self.N*self.Q), dtype=complex)
    
    def inv_T_s(self):
        temp = np.zeros((self.Q, self.N), dtype=complex)
        for q in range(self.Q):
            for n in range(self.N):
                temp[q,n] = self.ez[q,n] * self.sx[q,n] * self.sy[q,n]
        
        T_temp = np.copy(temp)

        row = np.array([i for i in range(self.N*self.Q)])
        column = np.array([i for i in range(self.N*self.Q)])
        
        value = np.array([0 for i in range(self.N*self.Q)], dtype=complex)
        for i in range(self.Q):
            value[i*self.N: (i+1)*self.N] = 1/T_temp[i, :]
        
        return csc_matrix((value, (row, column)), shape=(self.N*self.Q, self.N*self.Q), dtype=complex)

    def A_xx(self):
        row_main = np.array([i for i in range(self.N*self.Q)])
        column_main = np.array([i for i in range(self.N*self.Q)])
        value_main = np.array([], dtype=float)

        row_up = np.array([])
        column_up = np.array([])
        value_up =  np.array([], dtype=float)

        for q in range(self.Q):
            temp_row_up = np.array([i for i in range(q*self.N, (q+1)*self.N-1)])
            temp_column_up = np.array([i for i in range(q*self.N + 1, (q+1)*self.N)])

            temp_value_main = np.array([-self.t_xi[n]/self.d_X for n in range(self.N)])
            temp_value_up = np.array([self.t_xi[n]/self.d_X for n in range(1, self.N)])

            value_main = np.append(value_main, temp_value_main)
            value_up = np.append(value_up, temp_value_up)

            row_up = np.append(row_up, temp_row_up)
            column_up = np.append(column_up, temp_column_up)

        value_main[-1] = self.t_xi[self.N-1]/self.d_X
        row = np.append(row_main, row_up)
        column = np.append(column_main, column_up)
        value = np.append(value_main, value_up)


        return csc_matrix((value, (row, column)), shape=(self.N*self.Q, self.N*self.Q), dtype=complex)

    def A_yy(self):
        row_main = np.array([i for i in range(self.N*self.Q)])
        column_main = np.array([i for i in range(self.N*self.Q)])
        value_main = np.array([], dtype=float)

        row_up = np.array([i for i in range(self.N*(self.Q-1))])
        column_up = np.array([i for i in range(self.N, self.N*self.Q)])
        value_up = np.array([], dtype=float)

        for q in range(self.Q):
            temp_value_main = -self.t_xi/self.d_Y
            if q != self.Q-1:
                temp_value_up = self.t_xi/self.d_Y
                value_up = np.append(value_up, temp_value_up)

            value_main = np.append(value_main, temp_value_main)
            
        row = np.append(row_main, row_up)
        column = np.append(column_main, column_up)
        value = np.append(value_main, value_up)

        return csc_matrix((value, (row, column)), shape=(self.N*self.Q, self.N*self.Q), dtype=complex)
    
    def C_xx(self):
        row_main = np.array([i for i in range(self.N*self.Q)])
        column_main = np.array([i for i in range(self.N*self.Q)])
        value_main = np.array([], dtype=float)

        row_low = np.array([])
        column_low = np.array([])
        value_low = np.array([], dtype=float)

        for q in range(self.Q):
            temp_row_low = np.array([i for i in range(q*self.N + 1, (q+1)*self.N)])
            temp_column_low = np.array([i for i in range(q*self.N, (q+1)*self.N-1)])
            temp_value_main = np.array([self.t_xi[n]/self.d_X for n in range(self.N)])
            temp_value_low = np.array([-self.t_xi[n]/self.d_X for n in range(0, self.N-1)])

            value_main = np.append(value_main, temp_value_main)
            value_low = np.append(value_low, temp_value_low)

            row_low = np.append(row_low, temp_row_low)
            column_low = np.append(column_low, temp_column_low)
        
        row = np.append(row_main, row_low)
        column = np.append(column_main, column_low)
        value = np.append(value_main, value_low)

        return csc_matrix((value, (row, column)), shape=(self.N*self.Q, self.N*self.Q), dtype=complex)

    def C_yy(self):
        row_main = np.array([i for i in range(self.N*self.Q)])
        column_main = np.array([i for i in range(self.N*self.Q)])
        value_main = np.array([], dtype=float)

        row_low = np.array([i for i in range(self.N, self.N*self.Q)])
        column_low = np.array([i for i in range(self.N*(self.Q-1))])
        value_low = np.array([], dtype=float)

        for q in range(self.Q):
            if q < self.Q-1:
                temp_value_main = self.t_xi/self.d_Y
                temp_value_low = -self.t_xi/self.d_Y

                value_low = np.append(value_low, temp_value_low)
            else:
                temp_value_main = self.t_xi/self.d_Y

            value_main = np.append(value_main, temp_value_main)

        row = np.append(row_main, row_low)
        column = np.append(column_main, column_low)
        value = np.append(value_main, value_low)

        return csc_matrix((value, (row, column)), shape=(self.N*self.Q, self.N*self.Q), dtype=complex)
    
    @staticmethod
    def block_mat(A, B, C, D):
        _size_ = np.shape(A)[0]

        row_ul, col_ul = A.nonzero()
        val_ul = np.array([A[i, j] for i,j in zip(*A.nonzero())], dtype=complex)

        row_ur, col_ur = B.nonzero()
        val_ur = np.array([B[i, j] for i,j in zip(*B.nonzero())], dtype=complex)

        row_dl, col_dl = C.nonzero()
        val_dl = np.array([C[i, j] for i,j in zip(*C.nonzero())], dtype=complex)

        row_dr, col_dr = D.nonzero()
        val_dr = np.array([D[i, j] for i,j in zip(*D.nonzero())], dtype=complex)
        
        #changing the numeration for each element
        col_ur = col_ur + _size_
        row_dl = row_dl + _size_
        col_dr = col_dr + _size_
        row_dr = row_dr + _size_

        column = np.append(col_ul, np.append(col_ur, np.append(col_dl, col_dr)))
        row = np.append(row_ul, np.append(row_ur, np.append(row_dl, row_dr)))
        value = np.append(val_ul, np.append(val_ur, np.append(val_dl, val_dr)))

        return csc_matrix((value, (row, column)), shape=(2*_size_, 2*_size_), dtype=complex)

    def FDE(self, num: __init__):
        r'''Finite-difference eigensolver, which calculates eigenmodes and eigenvalues of curved waveguide
        
        Parameters
        --------
        ``num`` : int 
                number of calculating modes
        '''
        U_x, U_y, V_x, V_y = self.U_xx(self.N, self.Q, self.d_X), self.U_yy(self.N, self.Q, self.d_Y), self.V_xx(self.N, self.Q, self.d_X), self.V_yy(self.N, self.Q, self.d_Y)

        A_x = self.A_xx()
        A_y = self.A_yy()
        C_x = self.C_xx()
        C_y = self.C_yy()
        NPML = [150, 150, 60, 60]
        self.get_PML(NPML)
        S_xi, S_eta, inv_S_s = self.S_xi(self.N, self.Q, self.sx, self.sy, self.t_xi), self.S_eta(self.N, self.Q, self.sx, self.sy, self.t_xi), self.inv_S_s(self.N, self.Q, self.sx, self.sy)
        T_xi, T_eta, inv_T_s = self.T_xi(), self.T_eta(), self.inv_T_s()
            
        P_xx = -A_x @ inv_T_s @ V_y @ C_x @ inv_S_s @ U_y + (S_eta + A_x @ inv_T_s @ V_x) @ (T_xi + C_y @ inv_S_s @ U_y)
        P_xy = A_x @ inv_T_s @ V_y @ (T_eta + C_x @ inv_S_s @ U_x) - (S_eta + A_x @ inv_T_s @ V_x) @ C_y @ inv_S_s @ U_x
        P_yx = -(S_xi + A_y @ inv_T_s @ V_y) @ C_x @ inv_S_s @ U_y + A_y @ inv_T_s @ V_x @ (T_xi + C_y @ inv_S_s @ U_y)
        P_yy = (T_eta + C_x @ inv_S_s @ U_x) @ (S_xi + A_y @ inv_T_s @ V_y) - A_y @ inv_T_s @ V_x @ C_y @ inv_S_s @ U_x
        
        P = self.block_mat(P_xx, P_xy, P_yx, P_yy)

        ### EIGENVALUES ###
        vals, vecs = eigs(P, k=num, sigma=self.n_core**2)

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
        U_x, U_y, V_x, V_y = self.U_xx(self.N, self.Q, self.d_X), self.U_yy(self.N, self.Q, self.d_Y), self.V_xx(self.N, self.Q, self.d_X), self.V_yy(self.N, self.Q, self.d_Y)

        E_field = np.zeros((3*self.N*self.Q, mode_num), dtype=complex)

        C_x, C_y = self.C_xx(), self.C_yy()
        T_xi, T_eta, inv_T_s = self.T_xi(), self.T_eta(),  self.inv_T_s()
        inv_S_s = self.inv_S_s(self.N, self.Q, self.sx, self.sy)

        for i in range(1, mode_num + 1):
            h_x = (1/self.n_eff[i - 1]) * (C_x @ inv_S_s @ U_y @ self.vecs[:self.N*self.Q, i-1] - (T_eta + C_x @ inv_S_s @ U_x) @ self.vecs[self.N*self.Q : 2*self.N*self.Q, i-1])
            h_y = (1/self.n_eff[i - 1]) * ((T_xi + C_y @ inv_S_s @ U_y) @ self.vecs[:self.N*self.Q, i-1] - C_y @ inv_S_s @ U_x @ self.vecs[self.N*self.Q : 2*self.N*self.Q, i-1])

            e_z = 1j * inv_T_s @ (V_y @ h_x - V_x @ h_y)
            norm = np.sqrt((np.sum(np.abs(self.vecs[:self.N*self.Q, i-1])**2) + np.sum(np.abs(self.vecs[self.N*self.Q : 2*self.N*self.Q, i-1])**2) + np.sum(np.abs(e_z)**2)) * self.d_X * self.d_Y)
            
            vec_norm_x = self.vecs[:self.N*self.Q, i-1]/norm
            vec_norm_y = self.vecs[self.N*self.Q : 2*self.N*self.Q, i-1]/norm
            vec_norm_z = e_z/norm


            temp_vec = np.append(np.append(vec_norm_x, vec_norm_y), vec_norm_z, axis = 0)
            
            E_field[:, i-1] = temp_vec
        return E_field

    def draw_field(self, mode:int, scale:str):
        r'''Draw a density field projections for both electric and magnetic density profiles for orthonormal curvilinear axes ``xi``, ``eta`` and ``s`` in defined scale

        Parameters
        ---------
        ``mode`` : int 
                mode number 
        ``scale`` : str 
                scaling format for data visualization (`log` or `norm`)
        '''
        x_ticks = np.round(10**6 * np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N), 1)
        y_ticks = np.round(10**6 * np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q), 1)

        E_x, E_y, E_z, H_x, H_y, H_z = self.get_field(mode_num=mode)[:6]

        fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(24, 10),
                       gridspec_kw={
                           'width_ratios': [1, 1, 1],
                           'height_ratios': [5, 5],
                       'wspace': 0.2,
                       'hspace': 0.4})

        fig.suptitle(r'Field projections for ' + str(mode) +  r' mode : $n_{eff}$ =' + str(np.round(self.n_eff[mode-1], 3)) + r", $\kappa = $" + str(self.kappa) + r" $\mu m^{-1}$", fontsize=20)
        
        df_xi_E = pd.DataFrame(data=E_x, index=y_ticks, columns=x_ticks)
        df_eta_E = pd.DataFrame(data=E_y, index=y_ticks, columns=x_ticks)
        df_s_E = pd.DataFrame(data=E_z, index=y_ticks, columns=x_ticks)
        df_xi_H = pd.DataFrame(data=H_x, index=y_ticks, columns=x_ticks)
        df_eta_H = pd.DataFrame(data=H_y, index=y_ticks, columns=x_ticks)
        df_s_H = pd.DataFrame(data=H_z, index=y_ticks, columns=x_ticks)

        if scale == 'norm':
            axs[0,0].set_title(r'$E_{\xi}$', fontsize=15)
            sns.heatmap(np.abs(df_xi_E), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[0,0], cmap='jet')
            axs[0,0].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[0,0].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[0,0].set_yticklabels(y_ticks[::self.Q//7], rotation=0)

            
            axs[0,1].set_title(r'$E_{\eta}$', fontsize=15)
            sns.heatmap(np.abs(df_eta_E), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[0,1], cmap='jet')
            axs[0,1].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[0,1].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[0,1].set_yticklabels(y_ticks[::self.Q//7], rotation=0)
            
            axs[0,2].set_title(r'$E_{s}$', fontsize=15)
            sns.heatmap(np.abs(df_s_E), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[0,2], cmap='jet')
            axs[0,2].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[0,2].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[0,2].set_yticklabels(y_ticks[::self.Q//7], rotation=0)
            
            axs[1,0].set_title(r'$H_{\xi}$', fontsize=15)
            sns.heatmap(np.abs(df_xi_H), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[1,0], cmap='jet')
            axs[1,0].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[1,0].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[1,0].set_yticklabels(y_ticks[::self.Q//7], rotation=0)
            
            axs[1,1].set_title(r'$H_{\eta}$', fontsize=15)
            sns.heatmap(np.abs(df_eta_H), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[1,1], cmap='jet')
            axs[1,1].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[1,1].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[1,1].set_yticklabels(y_ticks[::self.Q//7], rotation=0)

            
            axs[1,2].set_title(r'$H_{s}$', fontsize=15)
            sns.heatmap(np.abs(df_s_H), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[1,2], cmap='jet')
            axs[1,2].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[1,2].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[1,2].set_yticklabels(y_ticks[::self.Q//7], rotation=0)

        
        elif scale == 'log':
            axs[0,0].set_title(r'$E_{\xi}$', fontsize=15)
            sns.heatmap(np.abs(df_xi_E), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[0,0], cmap='jet', norm=LogNorm())
            axs[0,0].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[0,0].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[0,0].set_yticklabels(y_ticks[::self.Q//7], rotation=0)

            
            axs[0,1].set_title(r'$E_{\eta}$', fontsize=15)
            sns.heatmap(np.abs(df_eta_E), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[0,1], cmap='jet', norm=LogNorm())
            axs[0,1].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[0,1].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[0,1].set_yticklabels(y_ticks[::self.Q//7], rotation=0)
            
            axs[0,2].set_title(r'$E_{s}$', fontsize=15)
            sns.heatmap(np.abs(df_s_E), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[0,2], cmap='jet', norm=LogNorm())
            axs[0,2].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[0,2].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[0,2].set_yticklabels(y_ticks[::self.Q//7], rotation=0)
            
            axs[1,0].set_title(r'$H_{\xi}$', fontsize=15)
            sns.heatmap(np.abs(df_xi_H), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[1,0], cmap='jet', norm=LogNorm())
            axs[1,0].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[1,0].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[1,0].set_yticklabels(y_ticks[::self.Q//7], rotation=0)
            
            axs[1,1].set_title(r'$H_{\eta}$', fontsize=15)
            sns.heatmap(np.abs(df_eta_H), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[1,1], cmap='jet', norm=LogNorm())
            axs[1,1].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[1,1].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[1,1].set_yticklabels(y_ticks[::self.Q//7], rotation=0)

            
            axs[1,2].set_title(r'$H_{s}$', fontsize=15)
            sns.heatmap(np.abs(df_s_H), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[1,2], cmap='jet', norm=LogNorm())
            axs[1,2].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[1,2].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[1,2].set_yticklabels(y_ticks[::self.Q//7], rotation=0)

        
        
        plt.show()

    def draw_E_field(self, mode:int, scale:str):
        r'''Draw a electric density field projections for orthonormal curvilinear axes ``xi``, ``eta`` and ``s`` in defined scale

        Parameters
        -----------
        ``mode`` : int
                mode number 
        ``scale`` :  str 
                scaling format for data visualization (`log` or `norm`)
        '''
        x_ticks = np.round(10**6 * np.linspace(-self.W/2 - self.delta_l, self.W/2 + self.delta_r, self.N), 1)
        y_ticks = np.round(10**6 * np.linspace(-self.H/2 - self.delta_d, self.H/2 + self.delta_u, self.Q), 1)

        E= self.get_E_field(mode_num=mode)
        E_x = np.zeros((self.Q, self.N), dtype=complex)
        E_y = np.zeros((self.Q, self.N), dtype=complex)
        E_z = np.zeros((self.Q, self.N), dtype=complex)
        for i in range(self.Q):
            E_x[i, :] = E[i*self.N:(i+1)*self.N, mode-1]
            E_y[i, :] = E[self.N*self.Q + i*self.N : self.N*self.Q + (i+1)*self.N, mode-1]
            E_z[i, :] = E[2*self.N*self.Q + i*self.N : 2*self.N*self.Q + (i+1)*self.N, mode-1]


        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(24, 10),
                       gridspec_kw={
                           'width_ratios': [1, 1, 1],
                           'height_ratios': [5],
                       'wspace': 0.2,
                       'hspace': 0.4})

        fig.suptitle(r'Field projections for ' + str(mode) +  r' mode : $n_{eff}$ =' + str(np.round(self.n_eff[mode-1], 3)) + r", $\kappa = $" + str(self.kappa) + r" $\mu m^{-1}$", fontsize=20)
        

        df_xi_E = pd.DataFrame(data=E_x, index=y_ticks, columns=x_ticks)
        df_eta_E = pd.DataFrame(data=E_y, index=y_ticks, columns=x_ticks)
        df_s_E = pd.DataFrame(data=E_z, index=y_ticks, columns=x_ticks)

        if scale == 'norm':
            axs[0].set_title(r'$E_{\xi}$', fontsize=15)
            sns.heatmap(np.abs(df_xi_E)**2, xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[0], cmap='jet')
            axs[0].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[0].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[0].set_yticklabels(y_ticks[::self.Q//7], rotation=0)

            
            axs[1].set_title(r'$E_{\eta}$', fontsize=15)
            sns.heatmap(np.abs(df_eta_E)**2, xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[1], cmap='jet')
            axs[1].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[1].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[1].set_yticklabels(y_ticks[::self.Q//7], rotation=0)
            
            axs[2].set_title(r'$E_{s}$', fontsize=15)
            sns.heatmap(np.abs(df_s_E)**2, xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[2], cmap='jet')
            axs[2].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[2].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[2].set_yticklabels(y_ticks[::self.Q//7], rotation=0)

        elif scale == 'log':
            axs[0].set_title(r'$E_{\xi}$', fontsize=15)
            sns.heatmap(np.abs(df_xi_E), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[0], cmap='jet', norm=LogNorm())
            axs[0].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[0].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[0].set_yticklabels(y_ticks[::self.Q//7], rotation=0)

            
            axs[1].set_title(r'$E_{\eta}$', fontsize=15)
            sns.heatmap(np.abs(df_eta_E), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[1], cmap='jet', norm=LogNorm())
            axs[1].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[1].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[1].set_yticklabels(y_ticks[::self.Q//7], rotation=0)
            
            axs[2].set_title(r'$E_{s}$', fontsize=15)
            sns.heatmap(np.abs(df_s_E), xticklabels=self.N//4, yticklabels=self.Q//7, ax=axs[2], cmap='jet', norm=LogNorm())
            axs[2].set_xlabel(r'$\xi$ coordinate, $\mu m$')
            axs[2].set_ylabel(r'$\eta$ coordinate, $\mu m$')
            axs[2].set_yticklabels(y_ticks[::self.Q//7], rotation=0)
        
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