import numpy as np                          # Para operações com vetores e matrizes
import matplotlib.pyplot as plt             # Para plotar gráficos
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg  # Para integrar o matplotlib com Tkinter
import scipy.integrate as integrate         # Para resolver integrais
import scipy.optimize as optimize           # Para encontrar mínimos, máximos, raízes e ajustar curvas
import tkinter as tk                        # Para criar interfaces gráficas com Tkinter
from tkinter import ttk                     # Para importar widgets de Tkinter (como o slider)

#Parâmetros
h = 6.62607015e-34  #Constate de Planck [J*s]
c = 3e8             #Velocidade da luz [m/s]
k = 1.380649e-23    #Constante de Boltzmann [J/K]

def transcendental_equation(x: float) -> float:
    """
    Função que representa a equação transcendental:
        f(x)  = 5*(e^x - 1) - x*e^x
    """
    return 5 * (np.expm1(x)) - x * np.exp(x)

def solve_transcendental_equation(initial_guess: float = 5.0, precision: int = 3) -> float:
    """
    Resolve a equação transcendental usando o método de Newton

    Parâmetros
    ----------
        initial_guess : float
            Chute inicial para o método de Newton. valor padrão = 5.0.
        precision : int 
            Número de casas decimais no resultado. valor padrão = 3.

    Retorna
    ----------
    float
        raiz aproximada da equação
    """
    root = optimize.newton(transcendental_equation, initial_guess)
    return np.round(root, precision)

x = solve_transcendental_equation()
b = h * c / (k * x) #Constante de Wien

def Lambda_Max(T: float) -> float:
    """
    Funçao que representa a lei de Wien:
        lambda_max = (hc)/kTx
    """
    return b / T

def dW_dlambda(wavelength: float, T: float) -> float:
    """
    Calcula a densidade espectral de potência (Lei de Planck) por comprimento de onda.

    Parâmetros
    ----------
    wavelength : float
        Comprimento de onda em metros.
    T : float
        Temperatura absoluta em Kelvin.

    Retorna
    -------
    float
        Potência espectral por unidade de comprimento de onda.
    """
    exponent = h * c / (wavelength * k * T)
    numerator = 2 * np.pi * h * c**2
    denominator = (wavelength**5) * (np.exp(exponent) - 1)
    
    return numerator / denominator

#Função universal x^3/(e^x - 1) usada na integração da Lei de Planck.
f = lambda x: x**3 / (np.expm1(x))

def W(T: float, wavelength = None) -> float:
    """
    Calcula a potência espectral integrada de um corpo negro entre dois comprimentos de onda.
    
    Parâmetros
    ----------
    T : float
        Temperatura absoluta (K).
    wavelength : list or tuple, optional
        [lamda_min, lambda_max] em metros. Padrão = [0, inf], representado como None.
        
    Retorno
    -------
    float
        Potência irradiada por unidade de área integrada no intervalo especificado [W/m²].
    """
    lambda_min, lambda_max = (0, np.inf) if wavelength is None else wavelength
    
    xmin = 0 if np.isinf(lambda_max) else h * c / (k * T * lambda_max)
    xmax = np.inf if lambda_min == 0 else h * c / (k * T * lambda_min)
    
    C = 2 * np.pi * (k**4) * (T**4) / ((c**2) * (h**3))
    int_F = integrate.quad(f, xmin, xmax)[0]
    
    return C * int_F

def W_total(T: float) -> float:
    """Potência total do corpo negro em T."""    
    return np.round(W(T, [1e-8, np.inf]), 0)

def wavelength_to_rgb(w: float) -> tuple:
    """
    Converte um comprimento de onda (em nm) para a cor RGB correspondente no espectro visível.
    
    Parâmetros
    ----------
    w : float
        Comprimento de onda em nanômetros (nm).
        
    Retorno
    -------
    tuple
        Cor correspondente no formato RGB, com valores entre 0 e 1.
    """
    # Verificar se o comprimento de onda está fora do espectro visível
    if w < 380 or w > 750:
        return (0, 0, 0)  # Fora do espectro visível

    # Inicializar as variáveis de cor
    r, g, b = 0, 0, 0

    # Definir a cor com base no comprimento de onda
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

    # Ajuste para comprimentos de onda fora do espectro visível
    if w < 420:
        factor = 0.3 + 0.7 * (w - 380) / (420 - 380)
    elif w > 700:
        factor = 0.3 + 0.7 * (780 - w) / (780 - 700)
    else:
        factor = 1

    # Retornar a cor RGB ajustada
    return (r * factor, g * factor, b * factor)



def update_plot(T: float) -> None:
    """
    Atualiza o gráfico de radiação de corpo negro com base na temperatura fornecida.
    
    Parâmetros
    ----------
    T : float
        Temperatura do corpo negro em Kelvin (K).
        
    Retorno
    -------
    None
    """
    # Gerar o vetor de comprimentos de onda (em metros)
    x = np.linspace(100e-9, 2000e-9, 1000)
    
    # Limpar o eixo do gráfico
    ax.clear()
    
    # Calcular o comprimento de onda máximo usando a Lei de Wien
    lambda_max = Lambda_Max(T)
    
    # Converter o comprimento de onda para nm (nanômetros)
    wavelength_nm = x * 1e9
    
    # Calcular as intensidades de radiação para cada comprimento de onda
    intensities = np.array([dW_dlambda(l, T) for l in x])
    y_max = intensities.max()  # Intensidade máxima para normalização
    
    # Inicializar a imagem RGB para o gráfico
    rgb_image = np.zeros((400, len(x), 3))
    
    # Preencher a imagem RGB com base no comprimento de onda
    for i, wl in enumerate(wavelength_nm):
        rgb_image[:, i] = wavelength_to_rgb(wl)
    
    # Criar a máscara para a imagem com base na intensidade
    mask = np.zeros_like(rgb_image[:, :, 0])
    for i, y_curve in enumerate(intensities):
        y_pixel_limit = int((y_curve / y_max) * rgb_image.shape[0])
        mask[:y_pixel_limit, i] = 1
    
    # Aplicar a máscara à imagem RGB
    rgb_masked = np.zeros_like(rgb_image)
    for c in range(3):
        rgb_masked[:, :, c] = rgb_image[:, :, c] * mask + 1 * (1 - mask)
    
    # Plotar a imagem RGB no gráfico
    ax.imshow(
        rgb_masked,
        extent=[wavelength_nm[0], wavelength_nm[-1], 0, y_max],
        aspect='auto',
        origin='lower'
    )
    
    # Plotar a curva de intensidade de radiação (Lei de Planck)
    ax.plot(wavelength_nm, dW_dlambda(x, T),
            color='black', 
            label=f'Temperatura Absoluta: {np.round(T, 2)}K')
    
    # Plotar a linha do comprimento de onda de máxima intensidade
    ax.plot(
        (lambda_max * 1e9, lambda_max * 1e9),
        (0, dW_dlambda(lambda_max, T)),
        '--', color='yellow', lw=7,
        label=f'Intensidade máxima {np.round(lambda_max * 1e9, 2)} nm'
    )
   
    # Plotar a intensidade total
    ax.plot(0, 0, 'ro', label=f'Intensidade total: {W_total(T):.2e} W/m^2')
    
    # Plotar a curva de intensidades no eixo
    ax.plot(wavelength_nm, intensities, color='black', linewidth=2)
    
    # Configurar os limites do gráfico
    ax.set_xlim(wavelength_nm[0], wavelength_nm[-1])  # Limites do eixo X (em nm)
    
    ax.set_ylim(0, y_max * 1.1)  # Limites do eixo Y (intensidade)
    
    # Adicionar a legenda
    ax.legend(loc=1)
    
    # Adicionar rótulos aos eixos
    ax.set_xlabel('Comprimento de onda [nm]', size=15)
    ax.set_ylabel(r'Intensidade [$\frac{W}{m^2}$]', size=15)
    
    # Adicionar a grade
    ax.grid()
    
    # Atualizar o canvas
    canvas.draw()

# Função para criar a interface com Tkinter
def create_gui() -> None:
    """
    Cria a interface gráfica utilizando o Tkinter para ajustar a temperatura e visualizar o espectro de radiação.
    """
    # Criar a janela principal do Tkinter
    root = tk.Tk()
    root.title("Espectro de radiação de corpo negro")
    
    root.columnconfigure(0, weight=1) 
    
    # Criar o slider para ajustar a temperatura
    slider = ttk.Scale(root, from_=200, to=11000, orient='horizontal', command=lambda val: update_plot(float(val)), length=700)
    slider.set(5800)  # Valor inicial de temperatura
    slider.grid(row=1, column=0, padx=15, pady=(0, 15), sticky='ew')

    # Criar o gráfico com Matplotlib
    global fig, ax, canvas
    fig, ax = plt.subplots(figsize=(10, 6), facecolor = 'lightgray')
    
    # Colocar o gráfico no canvas do Tkinter
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=(0, 10), sticky='nsew')
    
    root.rowconfigure(0, weight=1)

    # Inicializar o gráfico com a temperatura inicial
    update_plot(5800)

    # Rodar a interface Tkinter
    root.mainloop()


# Chama a função para criar a interface
create_gui()