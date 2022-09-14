# -*- coding: utf-8 -*-
"""
Gerador de superfície para uma Tabela de Basquete que sempre acerta na cesta

Arremessa a bola com uma posição e velocidade iniciais definidas

@author: Pena
"""
from funcoes import *

######## Configurações ########
visao = 1 # 1 = aberta, 2 = fechada, 3 = Lado Direito da Tabela
mostraArremesso = True
mostraCompleto = True
mostraRebatida = True
mostraNormal = True
mostraQuadra = True

"""
Arremesso Referência
PX0 = 4600
PY0 = 0
PZ0 = 1850

VX0 = -3941
VY0 = 0
VZ0 = 6851
"""

PX0 = 4600
PY0 = 0
PZ0 = 1850

VX0 = -3941 - 30
VY0 = 0 + 400
VZ0 = 6851 + 40

preparaPlot(visao, mostraQuadra)

####### Configura a Superfície da Tabela #######
# Carrega de um Arquivo
Dots = np.load('data/Dots2.npy')
Angles = np.load('data/Angles2.npy')
plotTabela(Dots, Angles)

Theta = Angles[:,:,0]
Phi = Angles[:,:,1]
maxDotsX = np.max(Dots[:,:,0])

######## Faz Arremesso ########
P0 = np.array([PX0, PY0, PZ0])
V0 = np.array([VX0, VY0, VZ0])
#print("\n\n deltaVZ0 {}, deltaVY0 {}".format(deltaVZ0, deltaVY0))

######## Faz arremesso #########
hit, grid, PF, VF = shoot(P0, V0, shootTimeRange, rebatida=False, 
                          Dots=Dots, Angles=Angles, maxDotsX=maxDotsX,
                          grafico=mostraArremesso, completo=mostraArremesso and mostraCompleto)
if hit: # Se acertou, reflete o arremesso e calcula o erro.
    i,j = grid
    theta, phi = Angles[i][j]
    VR0 = refletor(PF, VF, theta, phi, grafico=mostraNormal)
    PR0 = PF
    erro = shoot(PF, refletor(PF, VF, theta, phi, grafico=mostraRebatida and mostraNormal),
                  rebatidaTimeRange, rebatida=True, grafico=mostraRebatida)
    #print("Theta = {}, Phi = {}, Erro = {}".format(theta, phi, erro))
    plotSquare(i, j, Dots, projeta=False)
else:
    pass
    print("Arremesso PASSOU direto")

plt.show()
