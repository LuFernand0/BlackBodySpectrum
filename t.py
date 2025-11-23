import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import scipy.integrate as integrate 
import scipy.optimize as optimize 

# --- Constantes Físicas ---
# Parâmetros
H = 6.62607015e-34  # Constante de Planck [J*s]
C = 3e8             # Velocidade da luz [m/s]
K_B = 1.380649e-23  # Constante de Boltzmann [J/K]

# --- Solução para a Lei de Wien (Constante 'b') ---

def transcendental_equation(x: float) -> float:
    """
    Função que representa a equação transcendental para a Lei de Wien:
        f(x) = 5*(e^x - 1) - x*e^x
    """
    return 5 * (np.expm1(x)) - x * np.exp(x)

def solve_transcendental_equation(initial_guess: float = 5.0, precision: int = 4) -> float:
    """
    Resolve a equação transcendental usando o método de Newton para encontrar o termo hc/(\lambda_max k T).
    O valor da raiz (x) é aproximadamente 4.965.
    """
    try:
        # Usa o método de Newton para encontrar a raiz
        root = optimize.newton(transcendental_equation, initial_guess)
        return np.round(root, precision)
    except RuntimeError:
        # Caso o método não convirja, retorna um valor padrão
        return 4.9651 

X_ROOT = solve_transcendental_equation()
# Constante de Wien: b = lambda_max * T = h * c / (k * x)
WIEN_CONSTANT = H * C / (K_B * X_ROOT) 

# --- Funções Físicas ---

def Lambda_Max(T: float) -> float:
    """
    Função que representa a lei de Wien: lambda_max = b / T.
    """
    return WIEN_CONSTANT / T

def dW_dlambda(wavelength: float, T: float) -> float:
    """
    Calcula a densidade espectral de potência (Lei de Planck) por comprimento de onda.
    Unidade: W / (m^3) ou W / (m^2 * m)
    """
    if wavelength <= 0:
        return 0.0
        
    exponent = H * C / (wavelength * K_B * T)
    # 2 * pi * h * c^2: Usado para Radiância Espectral (Power per unit area, per unit solid angle, per unit wavelength)
    # ou Excitância Espectral (Power per unit area, per unit wavelength, 2 * pi * h * c^2)
    # Usaremos 2 * pi * h * c^2 para Excitância Espectral (W/m^3)
    numerator = 2 * np.pi * H * C**2
    denominator = (wavelength**5) * (np.exp(exponent) - 1)
    
    return numerator / denominator

# Função universal x^3/(e^x - 1) usada na integração da Lei de Planck.
F = lambda x: x**3 / (np.expm1(x))

def W(T: float, wavelength_range: tuple | None = None) -> float:
    """
    Calcula a potência espectral integrada de um corpo negro no intervalo.
    """
    lambda_min, lambda_max = (0, np.inf) if wavelength_range is None else wavelength_range
    
    # Mudança de variável para a integração
    xmin = 0 if np.isinf(lambda_max) else H * C / (K_B * T * lambda_max)
    xmax = np.inf if lambda_min == 0 else H * C / (K_B * T * lambda_min)
    
    C_integral = 2 * np.pi * (K_B**4) * (T**4) / ((C**2) * (H**3))
    
    # Tenta calcular a integral
    try:
        int_F = integrate.quad(F, xmin, xmax)[0]
    except (integrate.IntegrationWarning, ValueError):
        return 0.0
    
    return C_integral * int_F

def W_total(T: float) -> float:
    """Potência total do corpo negro em T."""
    # Integra de 10 nm até infinito para evitar singularidade em 0 e capturar quase toda a potência
    return np.round(W(T, [10e-9, np.inf]), 0)

def wavelength_to_rgb(w: float) -> tuple:
    """
    Converte um comprimento de onda (em nm) para a cor RGB correspondente no espectro visível.
    (Lógica customizada mantida do código original do usuário)
    """
    if w < 380 or w > 750:
        return (0, 0, 0)

    r, g, b = 0, 0, 0

    if 380 <= w < 440:
        r = -(w - 440) / (440 - 380)
        g = 0
        b = 1
    elif 440 <= w < 490:
        r = 0
        g = (w - 440) / (490 - 440)
        b = 1
    elif 490 <= w < 510:
        r = 0
        g = 1
        b = -(w - 510) / (510 - 490)
    elif 510 <= w < 580:
        r = (w - 510) / (580 - 510)
        g = 1
        b = 0
    elif 580 <= w < 645:
        r = 1
        g = -(w - 645) / (645 - 580)
        b = 0
    elif 645 <= w < 750:
        r = 1
        g = 0
        b = 0

    if w < 420:
        factor = 0.3 + 0.7 * (w - 380) / (420 - 380)
    elif w > 700:
        factor = 0.3 + 0.7 * (780 - w) / (780 - 700)
    else:
        factor = 1

    return (r * factor, g * factor, b * factor)

# --- Configurações de UI/Plot ---
T_MIN = 200.0   # Mínimo em Kelvin
T_MAX = 11000.0 # Máximo em Kelvin
T_INIT = 5800.0 # Temperatura inicial (aproximada do Sol)

# Vetor de comprimentos de onda para plotagem (de 100 nm a 2000 nm)
WAVELENGTH_M = np.linspace(100e-9, 2000e-9, 1000)
WAVELENGTH_NM = WAVELENGTH_M * 1e9

class BlackbodyApp:
    """
    Classe principal para gerenciar a interface gráfica e o estado do Matplotlib.
    """
    def __init__(self, master: tk.Tk):
        self.master = master
        master.title("Espectro de Radiação de Corpo Negro")
        
        # --- Variáveis de Estado ---
        self.current_T = T_INIT
        self.max_intensity = 0.0
        
        # --- Artistas do Matplotlib para Atualização Eficiente ---
        self.line_spectrum = None
        self.line_peak = None
        self.imshow_color = None
        
        self.setup_layout()
        self.setup_plot()
        self._update_plot_artists(self.current_T) # Inicialização

    def setup_layout(self):
        """Configura o layout da grade e os widgets Tkinter."""
        self.master.columnconfigure(0, weight=1)
        self.master.rowconfigure(0, weight=1)

        # 1. Gráfico (Matplotlib)
        fig, self.ax = plt.subplots(figsize=(10, 6), facecolor='lightgray')
        self.fig = fig
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas_widget = self.canvas.get_tk_widget()
        # Coloca na coluna 0, linha 0 e permite expansão
        self.canvas_widget.grid(row=0, column=0, padx=10, pady=(10, 10), sticky='nsew')
        
        # 2. Frame para o Slider e Rótulo (Lateral)
        control_frame = ttk.Frame(self.master)
        control_frame.grid(row=0, column=1, padx=10, pady=10, sticky='ns')
        
        # 3. Rótulo de Temperatura
        self.temp_label = tk.Label(control_frame, text=f'Temperatura: {T_INIT:.2f} K', font=('Arial', 12, 'bold'))
        self.temp_label.pack(pady=10)
        
        # 4. Slider Vertical (ttk.Scale)
        self.slider = ttk.Scale(
            control_frame, 
            from_=T_MAX, # Invertido para que o topo seja mais quente
            to=T_MIN,
            orient='vertical', 
            command=self.on_slider_move,
            length=400 
        )
        self.slider.set(T_INIT)
        self.slider.pack(pady=10)
        
        # 5. Rótulo de Potência Total
        self.power_label = tk.Label(control_frame, text=f'Potência Total: {W_total(T_INIT):.0f} W/m²', font=('Arial', 10))
        self.power_label.pack(pady=10)

    def setup_plot(self):
        """Inicializa as configurações estáticas do gráfico Matplotlib."""
        self.ax.set_xlabel('Comprimento de Onda [nm]', size=10)
        # Rótulo corrigido para W/m³ ou W/(m²·m)
        self.ax.set_ylabel(r'Densidade Espectral [$\frac{W}{m^3}$]', size=10)
        self.ax.grid(True, linestyle='--', alpha=0.6)
        self.ax.set_xlim(WAVELENGTH_NM[0], WAVELENGTH_NM[-1]) 
        self.ax.set_title('Espectro de Radiação de Corpo Negro', fontsize=14)
        
        # Cria placeholders para os artistas para serem atualizados depois
        self.line_spectrum, = self.ax.plot(WAVELENGTH_NM, np.zeros_like(WAVELENGTH_NM), color='black', linewidth=2, label='Curva de Planck')
        self.line_peak, = self.ax.plot([], [], '--', color='yellow', lw=5, label='λ max')
        
        # Inicializa o objeto imshow
        # O extent será corrigido na primeira atualização
        self.imshow_color = self.ax.imshow(
            np.zeros((10, len(WAVELENGTH_NM), 3)),
            extent=[WAVELENGTH_NM[0], WAVELENGTH_NM[-1], 0, 1],
            aspect='auto',
            origin='lower',
            interpolation='nearest'
        )
        
        self.ax.legend(loc='upper right')

    def on_slider_move(self, val: str):
        """Manipulador de eventos do slider."""
        T = float(val)
        self.current_T = max(T, T_MIN) # Garante que a temperatura não caia abaixo do mínimo
        self._update_plot_artists(self.current_T)

    def _update_plot_artists(self, T: float) -> None:
        """
        Atualiza os dados dos artistas Matplotlib existentes.
        Isso é muito mais eficiente do que ax.clear() e redesenhar.
        """
        # 1. Cálculo dos Dados
        intensities = np.array([dW_dlambda(l, T) for l in WAVELENGTH_M])
        y_max = intensities.max() * 1.1 # Novo limite Y

        # 2. Atualização dos Artistas
        
        # Curva de Planck
        self.line_spectrum.set_ydata(intensities)
        self.line_spectrum.set_label(f'T: {T:.2f} K')
        
        # Posição do Pico (Lei de Wien)
        lambda_max = Lambda_Max(T)
        peak_nm = lambda_max * 1e9
        peak_intensity = dW_dlambda(lambda_max, T)
        
        self.line_peak.set_data(
            (peak_nm, peak_nm),
            (0, peak_intensity)
        )
        self.line_peak.set_label(f'λ max: {peak_nm:.2f} nm')
        
        # Rótulos e Título
        self.temp_label.config(text=f'Temperatura: {T:.2f} K')
        self.power_label.config(text=f'Potência Total: {W_total(T):.0f} W/m²')
        self.ax.set_title(f'Espectro de Radiação de Corpo Negro (Pico em {peak_nm:.2f} nm)', fontsize=14)
        
        # 3. Atualização da Área Colorida (Espectro Visível)
        rgb_image = np.zeros((40, len(WAVELENGTH_NM), 3))
        for i, wl in enumerate(WAVELENGTH_NM):
            rgb_image[:, i] = wavelength_to_rgb(wl)
            
        # Criação da máscara de intensidade (mais simples do que a lógica original)
        mask = np.zeros_like(rgb_image[:, :, 0])
        for i, y_curve in enumerate(intensities):
            # Normaliza a intensidade para a altura da imagem (40 pixels)
            normalized_y = (y_curve / self.max_intensity) if self.max_intensity > 0 else 0 
            y_pixel_limit = int(normalized_y * rgb_image.shape[0] / 1.1)
            mask[:y_pixel_limit, i] = 1

        # A intensidade total pode mudar muito. Ajustamos a escala da curva de intensidade.
        if y_max > self.max_intensity:
            self.max_intensity = y_max
            
        # Aplica a imagem e ajusta a escala Y
        rgb_masked = rgb_image * mask[..., np.newaxis] 
        self.imshow_color.set_data(rgb_masked)
        self.imshow_color.set_extent([WAVELENGTH_NM[0], WAVELENGTH_NM[-1], 0, self.max_intensity])
        self.ax.set_ylim(0, self.max_intensity)
        
        # Redesenha a legenda e o canvas
        self.ax.legend(loc='upper right')
        self.canvas.draw_idle()

if __name__ == '__main__':
    root = tk.Tk()
    app = BlackbodyApp(root)
    root.mainloop()