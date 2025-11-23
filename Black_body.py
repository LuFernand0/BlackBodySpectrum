import numpy as np                  #Para operações com vetores e matrizes
import matplotlib.pyplot as plt     #Para plotar gráficos
import scipy.integrate as integrate #Para resolver integrais
import scipy.optimize as optimize   #Para encontrar minimos, máximos, raízes e ajusrar curvas

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