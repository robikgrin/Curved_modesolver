import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from mode_solver import rect_WG



def index_width_dep(str: rect_WG, num:int, widths):
    r'''Effective index dependency for all modes from the width of the core

    Parameters
    ------------
    ``num`` : int
            Number of modes for calculations
    ``widths`` : ndarray or list
            Values of core widths used for calculations

    Return
    ------
    ``n_eff`` : ndarray (num x len(widths))
            Values of effective indexes for all calculated modes

    '''
    n_eff = np.zeros((num, len(widths)), dtype=complex)
    W0 = keke.W
    for i in range(len(widths)):
        keke.W = keke.set_width(widths[i])
        keke.FDE(num)
        n_eff[:, i] = keke.n_eff
    keke.W = W0
    
    plt.figure(figsize=(10,10))
    for i in range(len(widths)):
        plt.scatter(widths*10**6, n_eff[i, :], linewidths=2)
        plt.plot(widths*10**6, n_eff[i, :], s = 20)
    plt.xlabel(r'Width $W$, $\mu m$')
    plt.ylabel(r'Effective index value $n_{eff}$')
    plt.title(r'Effective index from the width of core dependency ($n_{core} = $' + f'{keke.n_core}' + r', $n_{clad}$ = ' + f'{keke.n_clad}')
    plt.show()
    return n_eff