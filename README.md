# Espectro da Radiação de Corpo Negro

Este projeto é uma **simulaçao interativa** desenvolvida em Python para fins educacionais e científicos, realizado para a aula do de computação básica. O foco principal é a visualização das leis que regem a radiação térmica emitida por um corpo negro.

## Conceitos Físicos

A aplicação simula e demonstra as seguintes leis fundamentais da termodinâmica quântica:

- **Lei de Planck**: Calcula a intensidade da radiação emitida em diferentes comprimentos de onda ($\lambda$) para uma dada Temperatura ($T$).
- **Lei do Deslocamento de Wien**: Determina o comprimento de onda de máxima emissão ($\lambda_{\text{max}}$), que é inversamente proporcional à temperatura.
- **Lei de Stefan-Boltzmann**: Calcula a potência total irradiada (integral sob a curva de Planck).

## Tecnologias Utilizadas

O projeto utiliza o conjunto de bibliotecas do Python:
- `numpy`
- `scipy`
- `matplotlib`
- `tkinter`

## Como Rodar o Projeto

### Pré-requisitos
Certifique-se que você tem o Python instalado:
```bash
python --version
```
Crie um (ambiente virtual)[https://docs.python.org/pt-br/3/library/venv.html] e ative-o.

copie o repositório:
```Bash
git clone https://github.com/LuFernand0/BlackBodySpectrum.git
```
baixe as bibliotecas:
```Bash
pip install -r requirements.txt
```
(*As Bibliotecas `tkinter` e `tkk` geralmente vêm instaladas com o Python padrão.*)


### Executando o Script
Execute o Script no seu terminal:
```bash
python Black_body.py
```