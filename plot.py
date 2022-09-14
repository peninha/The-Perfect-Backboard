# -*- coding: utf-8 -*-
"""
Gerador de superfície para uma Tabela de Basquete que sempre acerta na cesta

Plota os pontos da tabela e suas normais

@author: Pena
"""
from funcoes import *

######## Configurações ########
visao = 3 # 1 = aberta, 2 = fechada, 3 = Lado Direito da Tabela
mostraNormals = True
mostraTabela = True
mostraQuadrados = False
projetaPontos = False
scale = 15

ax = preparaPlot(visao)
ax.view_init(elev=0., azim=0.01)
####### Configura a Superfície da Tabela #######
#Carrega de um Arquivo
Dots = np.load('data/Dots2.npy')
Angles = np.load('data/Angles2.npy')

if plotTabela:
    plotTabela(Dots, Angles, setas=mostraNormals, scale=scale)
for i in range(gridSizeX):
    for j in range(gridSizeY):
        if mostraQuadrados:
            plotSquare(i,j, Dots, projeta=projetaPontos) # Desenha os Quadrados Válidos


plt.show()
