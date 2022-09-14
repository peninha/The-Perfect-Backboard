# -*- coding: utf-8 -*-
"""
Spyder Editor

Programa por @Penadoxo
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d


######## Configurações ########
bolaRaio = 119
aroRaio = 450/2
aroHaste = 150
aroAltura = 3000
moscaRaio = 5
tabelaComprimento = 900
tabelaLargura = 600
focoAltura = 300
shotAltura = 1850
shotDistancia = 4600
deltaAngulo = 0.002
deltaT = 0.0005
shootTimeRange = 1.6
rebatidaTimeRange = 0.4
visao = 2 # 1 = aberta, 2 = fechada
g =  9807 #mm/s²


######## Funções ########

def colisor(P):
    global bolaRaio
    if P[0] <= bolaRaio:
        return [round(P[1]/10)+45, round((P[2]-3000)/10)]
    return False

def shoot(P0, V0, timeRange=2, rebatida=False, size=1):
    global deltaT, grid, PF, VF, aroP
    grid = 0
    X = []
    Y = []
    Z = []
    V = []
    for t in np.arange(0, timeRange, deltaT):
        X.append(P0[0] + V0[0]*t)
        Y.append(P0[1] + V0[1]*t)
        Z.append(P0[2] + V0[2]*t -g*t**2/2)
        
        if len(X) > 2:
            V.append((np.array([X[-1],Y[-1],Z[-1]]) -
                      np.array([X[-2],Y[-2],Z[-2]]))/deltaT) #calcula Vx, Vy e Vz
            if rebatida:
                if detector(Z[-1],V[-1][2]):
                    #print("tempo para cesta = {}".format(t))
                    ax.scatter3D(X, Y, Z, s=size)
                    return calcErro([X[-1],Y[-1],Z[-1]], aroP)
            else:
                test = colisor((X[-1], Y[-1], Z[-1]))
                if test: #### Colidiu ####
                    #print("tempo para bater na tabela = {}".format(t))
                    grid = test
                    VF = np.array(V[-1])
                    PF = np.array([X[-1], Y[-1], Z[-1]])
                    ax.scatter3D(X, Y, Z, s=size)
                    return True
    #ax.scatter3D(X, Y, Z, s=280)
    ax.scatter3D(X, Y, Z, s=size)
    return False #### Não Colidiu / Não Detectou ####

def refletor(P, V, theta, phi, seta = False):
    normal = calcNormal(theta, phi)
    if seta:    ### Desenha a Normal do bloco colidido ###
        ax.quiver(P[0], P[1], P[2], # <-- starting point of vector
                normal[0]*300, normal[1]*300, normal[2]*300, # <-- directions of vector
                color = 'red', alpha = .5, lw = 3)
    return (np.identity(3) - 2*np.multiply(normal, np.array([normal]).T)).dot(V)
    
def detector(z, Vz):
    global aroAltura
    if Vz <= 0:
        if z <= aroAltura:
            return True
    return False
    
def calcNormal(theta, phi):
    return np.array([np.cos(phi)*np.cos(theta), -np.cos(phi)*np.sin(theta),
                     np.sin(phi)])

def calcErro(P, Target):
    deltaX = Target[0] - P[0]
    deltaY = Target[1] - P[1]
    ax.quiver(P[0], P[1], P[2], # <-- starting point of vector
                deltaX, deltaY, 0, # <-- directions of vector
                color = 'black', alpha = .2, lw = 2)
    return np.sqrt(deltaX**2 + deltaY**2)
    
def Newton(erro):
    calcDerivada(erro)
    F0 = erro

def calcDerivada(erro):
    global deltaAngulo, grid
    F0 = erro
    theta, phi0, x = Grid[grid[0]][grid[1]]
    phi1 = phi0 + deltaAngulo
    F1 = testaAngulo(theta, phi1)
    
def testaAngulo(theta, phi, P0, Vmenos1):
    global testRange, deltaT
    pass

"""
    V0 = reflect(P0, Vmenos1, grid)
    
    for t in np.arange(0, testRange, deltaT):
        T = t - t0
        X.append(P0[0] + V0[0]*T)
        Y.append(P0[1] + V0[1]*T)
        Z.append(P0[2] + V0[2]*T -g*T**2/2)
        
        if len(X) > 2:
             V.append(calcV(np.array((X[-1], Y[-1], Z[-1])),
                            np.array((X[-2], Y[-2], Z[-2])),deltaT))
        
        test = colisor((X[-1], Y[-1], Z[-1]))
        if test and hit == False:
            grid = test
            V0 = reflect(np.array((X[-1], Y[-1], Z[-1])), V[-1], grid)
            hit = True
        
        if hit:
            if detector(Z[-1],V[-1][2]):
                erro = calcErro([X[-1],Y[-1],Z[-1]], aroP)
                while erro >= moscaRaio:
                    
                    print(erro)
                break
                #pass
"""

######## Desenha Tabela #########
aroP = np.array([aroHaste + aroRaio, 0, aroAltura])
aro2P = np.array([bolaRaio - (aroP[0]-bolaRaio), 0, aroAltura])
tabelaP = np.array([0, -tabelaComprimento/2, aroAltura])
tabela2P = np.array([bolaRaio, -tabelaComprimento/2, aroAltura])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
aro = plt.Circle(aroP, aroRaio, color='r', fill=0)
#aro2 = plt.Circle(aro2P, aroRaio, color='r', fill=0)
mosca = plt.Circle(aroP, moscaRaio, color='blue', fill=1)
tabela = plt.Rectangle((tabelaP[1],tabelaP[2]),
                       tabelaComprimento, tabelaLargura, color='b', fill=0)
tabela2 = plt.Rectangle((tabela2P[1],tabela2P[2]),
                       tabelaComprimento, tabelaLargura, color='b', fill=0)
ax.add_patch(aro)
ax.add_patch(mosca)
#ax.add_patch(aro2)
ax.add_patch(tabela)
ax.add_patch(tabela2)
art3d.pathpatch_2d_to_3d(aro, z=aroP[2], zdir="z")
art3d.pathpatch_2d_to_3d(mosca, z=aroP[2], zdir="z")
#art3d.pathpatch_2d_to_3d(aro2, z=aroP[2], zdir="z")
art3d.pathpatch_2d_to_3d(tabela, z=0, zdir="x")
art3d.pathpatch_2d_to_3d(tabela2, z=bolaRaio, zdir="x")

if visao == 1:    
    ax.set_xlim(0, 5000)
    ax.set_ylim(-2500, 2500)
    ax.set_zlim(0, 5000)
else:
    ax.set_xlim(-100, 850)
    ax.set_ylim(-450, 450)
    ax.set_zlim(2700, 3600)

######## Prepara o Grid ########
Grid = [[np.array([0,0,0]) for x in range(61)] for y in range(91)]
#Grid coordenadas = Theta (horizontal), Phi (vertical), X (profundidade)
for i in range (91):
    y = i*10 - 450
    x1 = aroP[0] - bolaRaio
    x2 = shotDistancia - bolaRaio
    theta = np.arctan2(y*(x1+x2),
                       np.sqrt(y**2+x2**2)*np.sqrt(y**2+x1**2)+x2*x1-y**2)
    #theta = 0
    for j in range (61):
        Grid [i][j] = np.array([theta, 0, 0])
        normal = bolaRaio*calcNormal(theta, 0)
        x = 0
        z = j*10 + 3000
        ######## Desenha as Normais do Grid ##########
        #ax.quiver(x, y, z, # <-- starting point of vector
            #normal[0], normal[1], normal[2], # <-- directions of vector
            #color = 'red', alpha = .5, lw = 1)


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

######## Modifica Arremesso de Referência #########
incrementoVx = 0
incrementoVy = 250
incrementoVz = 0
V0 = np.array([VX0+incrementoVx, VY0+incrementoVy, VZ0+incrementoVz])

######## Faz arremesso #########

hit = shoot(P0, V0, shootTimeRange, rebatida=False)
if hit: # Se acertou, reflete o arremesso e calcula o erro.
    theta, phi, z = Grid[grid[0]][grid[1]]
    VR0 = refletor(PF, VF, theta, phi, True)
    PR0 = PF
    erro = shoot(PR0, VR0, rebatidaTimeRange, rebatida=True)
    while erro >= moscaRaio:
        """
        VR = refletor(PF, VF, theta, phi+deltaAngulo)
        erroLinha = shoot(PR0, VR, rebatidaTimeRange, rebatida=True)
        derivada = (erroLinha-erro)/deltaAngulo
        phi2 = phi - erro/derivada
        print('Novo PHI = {}'.format(phi2))
        """
        phi2 = phi + deltaAngulo
        print('Novo PHI = {}'.format(phi2))
        VR2 = refletor(PF, VF, theta, phi2, True)
        erro2 = shoot(PR0, VR2, rebatidaTimeRange, rebatida=True)
        if erro2 >= erro:
            print("erro maior")
            break
        erro = erro2
        phi = phi2
        
        
else:
    print("Arremesso NÃO colidiu")


plt.show()


"""




X = []
Y = []
Z = []
V = []
t0 = 0
hit = False
for t in np.arange(0, timeRange, deltaT):
    T = t - t0
    X.append(P0[0] + V0[0]*T)
    Y.append(P0[1] + V0[1]*T)
    Z.append(P0[2] + V0[2]*T -g*T**2/2)
    
    if len(X) > 2:
         V.append(calcV(np.array((X[-1], Y[-1], Z[-1])),
                        np.array((X[-2], Y[-2], Z[-2])),deltaT))
    
    test = colisor((X[-1], Y[-1], Z[-1]))
    if test and hit == False:
        grid = test
        theta, phi, x = Grid[grid[0]][grid[1]]
        Vmenos1 = V[-1]
        V0 = reflect(np.array((X[-1], Y[-1], Z[-1])), V[-1], theta, phi, True)
        t0 = T
        hit = True
    
    if hit:
        if detector(Z[-1],V[-1][2]):
            erro = calcErro([X[-1],Y[-1],Z[-1]], aroP)
            #while erro >= moscaRaio:
                
            #    print(erro)
            break
            #pass        
    

#ax.scatter3D(X, Y, Z, s=280)
ax.scatter3D(X, Y, Z, s=1)

#ax.scatter3D(0, -450, 3600, s=280)
plt.show()
"""


