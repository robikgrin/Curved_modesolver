import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

def solve_1d_slab(d, n_core, n_clad, wl, polarization='TE'):
    """Решает трансцендентное уравнение для 1D планарного волновода."""
    k0 = 2 * np.pi / wl
    V = k0 * (d / 2) * np.sqrt(n_core**2 - n_clad**2)
    
    modes = []
    m = 0
    ratio = (n_core / n_clad)**2 if polarization == 'TM' else 1.0
    
    while True:
        # Строгое физическое условие существования моды m
        if V <= m * np.pi / 2:
            break
            
        u_min = m * np.pi / 2 + 1e-6
        u_max = min((m + 1) * np.pi / 2 - 1e-6, V - 1e-6)
        
        # Защита от схлопывания интервала
        if u_min >= u_max:
            break
            
        def eq(u):
            w = np.sqrt(V**2 - u**2)
            if m % 2 == 0: # Четные моды
                return u * np.tan(u) - ratio * w
            else:          # Нечетные моды
                return -u / np.tan(u) - ratio * w
                
        # ЗАЩИТА ОТ ЧИСЛЕННЫХ АРТЕФАКТОВ НА ГРАНИЦЕ ОТСЕЧКИ
        if np.sign(eq(u_min)) == np.sign(eq(u_max)):
            break
                
        # Теперь метод brentq сработает безупречно
        res = root_scalar(eq, bracket=[u_min, u_max], method='brentq')
        
        if res.converged:
            n_eff = np.sqrt(n_core**2 - (res.root / (k0 * d / 2))**2)
            modes.append(n_eff)
            
        m += 1
        
    return modes

def get_2d_TE_modes(W, H, n_core, n_clad, wl):
    """
    Метод Эффективного Показателя Преломления (EIM) для прямоугольного волновода.
    Ищем 2D TE-моды (основное электрическое поле Ex - по горизонтали).
    """
    # ШАГ 1: Решаем вертикальный срез. 
    # Поле Ex параллельно границам H, поэтому для вертикали это TE-задача.
    n_slabs = solve_1d_slab(H, n_core, n_clad, wl, polarization='TE')
    
    if len(n_slabs) == 0:
        return [] # Волновод слишком тонкий, не держит свет
        
    n_slab_fundamental = n_slabs[0] # Берем фундаментальную вертикальную моду
    
    # ШАГ 2: Решаем горизонтальный срез.
    # Поле Ex перпендикулярно границам W, поэтому для горизонтали это TM-задача.
    final_modes = solve_1d_slab(W, n_slab_fundamental, n_clad, wl, polarization='TM')
    
    return final_modes

# =====================================================================
# ИСПОЛЬЗОВАНИЕ
# =====================================================================
wl = 1.55
n_core = 3.4755
n_clad = 1.444
H = 0.22
W_test = 1.0 # Твои 2 микрона

modes = get_2d_TE_modes(W_test, H, n_core, n_clad, wl)

print(f"=====================================================")
print(f" Анализ волновода: W = {W_test} мкм, H = {H} мкм")
print(f"=====================================================")
print(f"Максимальное количество TE мод: {len(modes)}")
for i, n_eff in enumerate(modes):
    print(f"  TE_{i} мода: n_eff = {n_eff:.4f}")
print(f"=====================================================\n")

# --- БОНУС: Строим график зависимости числа мод от ширины W ---
W_array = np.linspace(0.4, 3.0, 300)
mode_counts = [len(get_2d_TE_modes(w, H, n_core, n_clad, wl)) for w in W_array]

plt.figure(figsize=(10, 5), dpi=150)
plt.plot(W_array, mode_counts, color='navy', linewidth=2)
plt.fill_between(W_array, mode_counts, alpha=0.2, color='royalblue')
plt.axvline(W_test, color='red', linestyle='--', label=f'W = {W_test} мкм (Мод: {len(modes)})')
plt.title('Максимальное число распространяемых TE-мод (EIM)', fontsize=14)
plt.xlabel('Ширина волновода W, мкм', fontsize=12)
plt.ylabel('Количество мод', fontsize=12)
plt.yticks(range(0, max(mode_counts)+2))
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.show()