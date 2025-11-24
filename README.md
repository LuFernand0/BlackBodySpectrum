# 🔬 Espectro da Radiação de Corpo Negro

Este projeto é uma simulação interativa desenvolvida em **Python** para fins educacionais e científicos.  
Foi realizado como um projeto prático para a disciplina de **Computação Básica**.

O objetivo é fornecer uma **ferramenta visual** para estudar as leis que regem a radiação térmica emitida por um corpo negro, mostrando como o espectro de luz e a potência total variam com a temperatura.

---

## 🌟 Visão Geral e Funcionalidades

A aplicação possui uma **interface gráfica** (tkinter) que exibe o gráfico do espectro (matplotlib) e permite o controle da temperatura através de um *slider*.

### Funções Principais

| Funcionalidade | Descrição |
|---------------|-----------|
| **Simulação Dinâmica** | A curva de Planck é redesenhada em tempo real conforme a temperatura é ajustada. |
| **Controle de Temperatura** | Slider para variar a Temperatura Absoluta (T) na faixa de **200 K a 11000 K**. |
| **Pico de Wien** | Marcação visual de `λ_max`, demonstrando a **Lei do Deslocamento de Wien**. |
| **Espectro de Cores** | A faixa visível (380–750 nm) é colorida, mostrando a mudança de cor do corpo negro conforme T aumenta. |
| **Potência Total** | Exibe a potência total irradiada (Lei de Stefan–Boltzmann) em notação científica. |

---

## 💡 Conceitos Físicos

A simulação representa três leis fundamentais da física da radiação térmica:

![Espectro]{\Graficos\grafico1.png}


### **Lei de Planck**
Calcula a intensidade da radiação emitida em diferentes comprimentos de onda (λ) para uma dada temperatura (T).

### **Lei do Deslocamento de Wien**
Define o comprimento de onda de máxima emissão:

$$
\lambda_{\text{max}} \propto \frac{1}{T}
$$

### **Lei de Stefan–Boltzmann**
A potência total irradiada por unidade de área é dada por:

$$
W \propto T^4
$$

---

## 🔧 Tecnologias Utilizadas

| Categoria | Biblioteca | Uso no Projeto |
|----------|------------|----------------|
| **Computação Científica** | numpy | Vetores, constantes físicas, funções matemáticas |
| **Análise Numérica** | scipy | Método de Newton e integrações |
| **Visualização** | matplotlib | Plot do espectro e interface gráfica |
| **Interface Gráfica** | tkinter / ttk | Janela e slider interativo |

---

## 🚀 Como Rodar o Projeto

### **Pré-requisitos**
Certifique-se de ter **Python 3.x** instalado.

---

### **1. Clonar o Repositório e Configurar o Ambiente**

```bash
# 1. Copiar o repositório
git clone https://github.com/LuFernand0/BlackBodySpectrum.git
cd BlackBodySpectrum

# 2. (Opcional) Criar e ativar um ambiente virtual
python -m venv venv

# Linux
source venv/bin/activate

# Windows
venv\Scripts\activate
```

### **2. Instalar Dependências**

```bash
pip install -r requirements.txt
```

### **3. Executar**
```bash
python Black_body.py
```
