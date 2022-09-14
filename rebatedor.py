# -*- coding: utf-8 -*-
"""
Gerador de superfície para uma Tabela de Basquete que sempre acerta na cesta

Algoritmo que rebate os lançamentos e calcula as normais ideais em cada ponto

@author: Pena
"""
from funcoes import *

######## Configurações ########
visao = 2 # 1 = aberta, 2 = fechada, 3 = Lado Direito da Tabela
mostraArremesso = True
mostraCompleto = False
mostraRebatida = True
mostraNormal = False
mostraCorrecao = False
mostraQuadra = True
#VXZ0range = [-141,124,3]
#VY0range = [0,420,5]
#VXZ0range = range(0,101,50)
VXZ0range = range(0,101,50)
#VY0range = range(,301,50)
VY0range = range(000,301,50)


ax = preparaPlot(visao, mostraQuadra)
####### Configura a Superfície da Tabela #######
#Carrega de um Arquivo
Dots = np.load('data/Dots0.npy')
Angles = np.load('data/Angles0.npy')

# Ou gera os pontos
"""
Dots = np.zeros((gridSizeX,gridSizeY,3))
Angles = np.zeros((gridSizeX,gridSizeY,2))
for i in range(gridSizeX):
    for j in range(gridSizeY):
        Dots[i][j] = grid2Coords(i, j)
        Angles[i][j] = np.array([0, 0])
        #plotNormal(i, j, Angles, x = Dots[i][j][0])
"""

#plotTabela(Dots, Angles, setas=True)

######## Configura Variáveis Iniciais #########
AnglesAmostra = np.zeros((gridSizeX,gridSizeY),)
AnglesAmostra = [[ [] for y in range(gridSizeY)] for x in range(gridSizeX)]
Theta = Angles[:,:,0]
Phi = Angles[:,:,1]
maxDotsX = np.max(Dots[:,:,0])

######## Calcula Arremesso de Referência #########
P0 = np.array([shotDistancia, 0, shotAltura])
P1 = np.array([bolaRaio, 0, aroAltura+focoAltura])
P2 = np.array([aro2P[0], 0, aroAltura])
DeltaX1 = P1[0]-P0[0]
DeltaX2 = P2[0]-P0[0]
DeltaZ1 = P1[2]-P0[2]
DeltaZ2 = P2[2]-P0[2]
VX0 = -np.sqrt(g*(DeltaX1 - DeltaX2)/(2*(DeltaZ2/DeltaX2 - DeltaZ1/DeltaX1)))
t1 = DeltaX1/VX0
t2 = t1*DeltaX2/DeltaX1
VY0 = 0
VZ0 = DeltaZ1/t1 + g*t1/2

######## Loop de Arremessos ########
for deltaVZ0 in VXZ0range:
    for deltaVY0 in VY0range:
        ######## Modifica Arremesso de Referência ########
        deltaVX0 = -deltaVZ0
        V0 = np.array([VX0, VY0, VZ0]) + [deltaVX0, deltaVY0, deltaVZ0]
        #print("\n\n deltaVZ0 {}, deltaVY0 {}".format(deltaVZ0, deltaVY0))
        
        ######## Faz arremesso #########
        hit, grid, PF, VF = shoot(P0, V0, shootTimeRange, rebatida=False, 
                                  Dots=Dots, Angles=Angles, maxDotsX=maxDotsX,
                                  grafico=mostraArremesso, completo=mostraArremesso and mostraCompleto)
        if hit: # Se acertou, reflete o arremesso e calcula o erro.
            i,j = grid
            theta, phi = Angles[i][j]
            VR0 = refletor(PF, VF, theta, phi, grafico=mostraNormal)
            #desvio = grid2Coords(i, j) - PF
            #PR0 = PF + desvio
            #PR0 = Dots[i][j] + bolaRaio*calcNormal(theta, phi)
            PR0 = PF
            erro = shoot(PR0, VR0, rebatidaTimeRange, rebatida=True, grafico=mostraRebatida)
            #print("Theta = {}, Phi = {}, Erro = {}".format(theta, phi, erro))
            thetaTeste = theta
            phiTeste = phi
            tries = 1
            while erro >= moscaRaio and tries < 10:
                ########### Corrige Phi e Theta até acertar na Mosca #############
                a = np.radians(-45) #-35 graus
                b = np.radians(45) #+35 graus
                phiTeste = gss(a, b, PF, VF, thetaTeste, phiTeste, var="phi", tol=1e-5)
                erro = shoot(PF, refletor(PF, VF, thetaTeste, phiTeste, grafico=mostraCorrecao and mostraNormal),
                              rebatidaTimeRange, rebatida=True, grafico=mostraCorrecao)
                #print("Novo Phi = {}, Erro = {}".format(phi2, erro))
                a = np.radians(0)
                b = np.radians(90) #45 graus
                thetaTeste = gss(a, b, PF, VF, thetaTeste, phiTeste, var="theta", tol=1e-5)
                erro = shoot(PF, refletor(PF, VF, thetaTeste, phiTeste, grafico=mostraCorrecao and mostraNormal),
                              rebatidaTimeRange, rebatida=True, grafico=mostraCorrecao)
                tries += 1
            if tries >= 10:
                print("Solução NÃO encontrada")
            else:
                Angles[i][j] = np.array([thetaTeste, phiTeste])
                AnglesAmostra[i][j].append(Angles[i][j])
                print("\nGrid = {}".format(grid))
                #print("Theta = {}, Phi = {}, Erro = {}".format(thetaTeste, phiTeste, erro))
        else:
            pass
            print("Arremesso PASSOU direto")

########## Faz a média dos valores, atualiza no Grid e plota Normais e Quadrados ###########
try:
    AnglesNovo = np.load("data/AnglesNovo.npy")
except FileNotFoundError:
    AnglesNovo = np.zeros((gridSizeX,gridSizeY,2))

for i in range(gridSizeX):
    for j in range(gridSizeY):
        if AnglesAmostra[i][j]:
            print ("grid {} {} = {}".format(i, j, len(AnglesAmostra[i][j])))
            AnglesNovo[i][j] = np.sum(AnglesAmostra[i][j], axis=0)/(len(AnglesAmostra[i][j]))
            #plotNormal(i, j, AnglesNovo, x=bolaRaio)  # Desenha as Normais do Grid #
            plotSquare(i,j, Dots, projeta=False) # Desenha os Quadrados Válidos

np.save('data/AnglesNovo', AnglesNovo)

plt.show()
