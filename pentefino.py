# -*- coding: utf-8 -*-
"""
Gerador de superfície para uma Tabela de Basquete que sempre acerta na cesta

Faz o pente fino nos pontos da tabela para ver quais ainda não foram sondados

@author: Pena
"""

from funcoes import *

visao = 3 # 1 = aberta, 2 = fechada, 3 = Lado Direito da Tabela
preparaPlot(visao)
Dots = np.load('data/Dots0.npy')
AnglesNovo = np.load("data/AnglesNovo.npy")

for i in range(gridSizeX):
    for j in range(gridSizeY-1,-1,-1):
        if AnglesNovo[i][j][0] == 0 and AnglesNovo[i][j][1] == 0 :
            print("grid {} {}".format(i, j))
            plotSquare(i,j, Dots, projeta=False) # Desenha os Quadrados faltantes