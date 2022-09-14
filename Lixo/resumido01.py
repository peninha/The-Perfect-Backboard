# -*- coding: utf-8 -*-
"""
Gerador de superfície para uma Tabela de Basquete que sempre acerta na cesta
Somente alguns algoritmos úteis

Programa por @Penadoxo
"""


import math
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from scipy import interpolate

######## Configurações ########
visao = 3 # 1 = aberta, 2 = fechada, 3 = Lado Direito da Tabela
#Arremessos = 500
VXZ0range = [-141,124,3]
VY0range = [0,420,5]
#VXZ0range = [109,124,3] #Varredura Inferior
#VY0range = [0,420,5]
#VXZ0range = [-141,-126,3] #Varredura Superior
#VY0range = [0,420,5]
#VXZ0range = [-140,140,3] #Varredura Centro
#VY0range = [0,20,10]
#VXZ0range = [-140,140,3] #Varredura Direita
#VY0range = [370,420,5]
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
SigmaVx = 80
SigmaVy = 200
SigmaVz = 80
deltaAngulo = 0.002
deltaT = 0.0005
shootTimeRange = 3
rebatidaTimeRange = 1
g =  9807 #mm/s²
gridSizeX = 91
gridSizeY = 61

SquareComprimento = tabelaComprimento/(gridSizeX-1)
SquareLargura = tabelaLargura/(gridSizeY-1)

######## Funções ########
def colisor(P):
    global bolaRaio, aroAltura
    if P[0] <= bolaRaio:
        return [round(P[1]/10)+45, round((P[2]-3000)/10)]
    return False

def shoot(P0, V0, timeRange=2, rebatida=False, grafico=False, size=1):
    global deltaT, PF, VF, aroP
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
                    if grafico:
                        ax.scatter3D(X, Y, Z, s=size)
                    return calcErro([X[-1],Y[-1],Z[-1]], aroP, seta=grafico)
            else:
                test = colisor((X[-1], Y[-1], Z[-1]))
                if test and test[0] in range(gridSizeX) and test[1] in range(gridSizeY): #### Colidiu ####
                    #print("tempo para bater na tabela = {}".format(t))
                    grid = test
                    VF = np.array(V[-1])
                    PF = np.array([X[-1], Y[-1], Z[-1]])
                    if grafico:
                        ax.scatter3D(X, Y, Z, s=size)
                    return True, grid
    if grafico:
        ax.scatter3D(X, Y, Z, s=size)
    return False, False #### Não Colidiu / Não Detectou ####

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

def calcErro(P, Target, seta=False):
    deltaX = Target[0] - P[0]
    deltaY = Target[1] - P[1]
    if seta:
        ax.quiver(P[0], P[1], P[2], # <-- starting point of vector
                deltaX, deltaY, 0, # <-- directions of vector
                color = 'black', alpha = .2, lw = 2)
    return np.sqrt(deltaX**2 + deltaY**2)

def grid2Coords(gridX, gridY):
    global  gridSizeX, gridSizeY, tabelaComprimento, tabelaLargura, focoAltura, aroAltura
    y = (gridX - (gridSizeX-1)/2) * tabelaComprimento/(gridSizeX-1)
    z = (gridY - (gridSizeY-1)/2) * tabelaLargura/(gridSizeY-1) + focoAltura + aroAltura
    x = 0
    return np.array([x, y, z])

def coords2Grid(P):
    global  gridSizeX, gridSizeY, tabelaComprimento, tabelaLargura, focoAltura, aroAltura
    x = P[1]
    y = P[2] - focoAltura - aroAltura
    gridX = (2*x + tabelaComprimento)*(gridSizeX-1) / (2*tabelaComprimento)
    gridY = (2*y + tabelaLargura)*(gridSizeY-1) / (2*tabelaLargura)
    return np.array([int(gridX), int(gridY)])

def grid2SquareCoords(gridX, gridY):
    global SquareComprimento, SquareLargura
    P = grid2Coords(gridX, gridY)
    x = P[1] - SquareComprimento/2
    y = P[2] - SquareLargura/2
    return np.array([x,y])

invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2
def gss(a, b, PF, VF, theta, phi, var="phi", tol=1e-5):
    """Golden-section search.
    Given a function f with a single local minimum in
    the interval [a,b], gss returns a subset interval
    [c,d] that contains the minimum with d-c <= tol.
    """
    (a, b) = (min(a, b), max(a, b))
    h = b - a
    if h <= tol:
        return (a, b)
    # Required steps to achieve tolerance
    n = int(math.ceil(math.log(tol / h) / math.log(invphi)))
    c = a + invphi2 * h
    d = a + invphi * h
    if var == "phi":
        fc = shoot(PF, refletor(PF, VF, theta, c, seta = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
        fd = shoot(PF, refletor(PF, VF, theta, d, seta = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
    else:
        fc = shoot(PF, refletor(PF, VF, c, phi, seta = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
        fd = shoot(PF, refletor(PF, VF, d, phi, seta = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
    yc = fc
    yd = fd
    for k in range(n-1):
        if yc < yd:  # yc > yd to find the maximum
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            if var == "phi":
                fc = shoot(PF, refletor(PF, VF, theta, c, seta = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
            else:
                fc = shoot(PF, refletor(PF, VF, c, phi, seta = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
            yc = fc
        else:
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            if var == "phi":
                fd = shoot(PF, refletor(PF, VF, theta, d, seta = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
            else:
                fd = shoot(PF, refletor(PF, VF, d, phi, seta = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
            yd = fd
    if yc < yd:
        #return (a, d)
        return (a + d)/2
    else:
        #return (c, b)
        return (c + b)/2

######### Configura o PyPlot ##########
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho', facecolor='#AAAAAA')
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
if visao == 1:    
    ax.set_xlim(0, 5000)
    ax.set_ylim(-2500, 2500)
    ax.set_zlim(0, 5000)
elif visao == 2:
    ax.set_xlim(-100, 800)
    ax.set_ylim(-450, 450)
    ax.set_zlim(2700, 3600)
else:
    ax.set_xlim(-100, 350)
    ax.set_ylim(0, 450)
    ax.set_zlim(3100, 3550)
ax.axis("off")

######## Desenha Tabela e Aro #########
aroP = np.array([aroHaste + aroRaio, 0, aroAltura])
aro2P = np.array([bolaRaio - (aroP[0]-bolaRaio), 0, aroAltura])
tabelaP = np.array([0, -tabelaComprimento/2, aroAltura])
tabela2P = np.array([bolaRaio, -tabelaComprimento/2, aroAltura])
aro = plt.Circle(aroP, aroRaio, color='r', fill=0)
mosca = plt.Circle(aroP, moscaRaio, color='blue', fill=1)
tabela = plt.Rectangle((tabelaP[1],tabelaP[2]),
                       tabelaComprimento, tabelaLargura, color='w', fill=1)
tabela2 = plt.Rectangle((tabela2P[1],tabela2P[2]),
                       tabelaComprimento, tabelaLargura, color='b', fill=0)
ax.add_patch(aro)
ax.add_patch(mosca)
ax.add_patch(tabela)
ax.add_patch(tabela2)
art3d.pathpatch_2d_to_3d(aro, z=aroP[2], zdir="z")
art3d.pathpatch_2d_to_3d(mosca, z=aroP[2], zdir="z")
art3d.pathpatch_2d_to_3d(tabela, z=0, zdir="x")
art3d.pathpatch_2d_to_3d(tabela2, z=bolaRaio, zdir="x")


###############################################
X = np.arange(45, 91, 1)
Y = np.arange(0, 61, 1)

delx = []
dely = []
for i in range(45,91):
    for j in range(2,61):
        if not isinstance(GridAmostras[i][j][0], np.ndarray):
            print ("grid {} {}".format(i, j))
            delx.append(i-45)
            dely.append(j-2)
            P = grid2SquareCoords(i, j)
            gridSquare = plt.Rectangle(P, SquareComprimento, SquareLargura, color='blue', fill=1)
            ax.add_patch(gridSquare)
            art3d.pathpatch_2d_to_3d(gridSquare, z=bolaRaio, zdir="x")
plt.show()
            
delx = list(dict.fromkeys(delx))
dely = list(dict.fromkeys(dely))
X = np.delete(X, delx)
Y = np.delete(Y, dely)
#xx, yy = np.meshgrid(x, y)

Theta = [[np.array([0]) for x in range(len(X))] for y in range(len(Y))]
Phi = [[np.array([0]) for x in range(len(X))] for y in range(len(Y))]
ThetaPhi = [[np.array([0,0,0]) for x in range(len(X))] for y in range(len(Y))]
for j in range(len(Y)):
    for i in range (len(X)):
        Theta[j][i] = Grid[X[i]][Y[j]][0]
        Phi[j][i] = Grid[X[i]][Y[j]][1]
        ThetaPhi[j][i] = Grid[X[i]][Y[j]]
        
fTheta = interpolate.interp2d(X, Y, Theta, kind='cubic')
fPhi = interpolate.interp2d(X, Y, Phi, kind='cubic')

def f(x, y):
    global fTheta, fPhi
    return np.concatenate((fTheta(x,y),fPhi(x,y)))




"""
########## Faz a média dos valores, atualiza no Grid e plota Normais e Quadrados ###########
for i in range(len(Grid)):
    for j in range(len(Grid[0])):
        if isinstance(GridAmostras[i][j][0], np.ndarray):
            #print ("grid {} {} = {}".format(i, j, len(GridAmostras[i][j])-1))
            Grid[i][j] = np.sum(GridAmostras[i][j], axis=0)/(len(GridAmostras[i][j])-1)
            
            ######## Desenha as Normais do Grid ##########
            P = grid2Coords(i, j)
            normal = bolaRaio/4*calcNormal(Grid[i][j][0],Grid[i][j][1])
            ax.quiver(P[0]+bolaRaio, P[1], P[2], # <-- starting point of vector
            normal[0], normal[1], normal[2], # <-- directions of vector
            color = 'red', alpha = 0.5, lw = 2)
            
            ######## Desenha os Quadrados #######
            P = grid2SquareCoords(i, j)
            gridSquare = plt.Rectangle(P, SquareComprimento, SquareLargura, color='blue', fill=1)
            ax.add_patch(gridSquare)
            art3d.pathpatch_2d_to_3d(gridSquare, z=bolaRaio, zdir="x")
            
plt.show()
"""