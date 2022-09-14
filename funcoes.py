# -*- coding: utf-8 -*-
"""
Gerador de superfície para uma Tabela de Basquete que sempre acerta na cesta

Funções gerais

@author: Pena
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Arc

######## Configurações ########
visao = 1 # 1 = aberta, 2 = fechada, 3 = Lado Direito da Tabela

bolaRaio = 119
aroRaio = 450/2
aroHaste = 150
#aroHaste = 110
aroAltura = 3000
moscaRaio = 5
tabelaComprimento = 900
tabelaLargura = 550
tabelaAltura = 3050
focoAltura = 300
shotAltura = 1850
shotDistancia = 4600
SigmaVx = 80
SigmaVy = 200
SigmaVz = 80
deltaAngulo = 0.002
#deltaT = 0.0001
deltaT = 0.0005
shootTimeRange = 3
rebatidaTimeRange = 1
g =  9807 #mm/s²
gridSizeX = 46
#gridSizeY = 61
gridSizeY = 56
quadraLargura = 15000
quadraComprimento = 11000
garrafaoLargura = 4850
garrafaoComprimento = 5800
garrafaoCirculoDiametro = 3600
mostraQuadra = False


squareComprimento = tabelaComprimento/(gridSizeX-1)/2
squareLargura = tabelaLargura/(gridSizeY-1)
aroP = np.array([aroHaste + aroRaio, 0, aroAltura])
aro2P = np.array([bolaRaio - (aroP[0]-bolaRaio), 0, aroAltura])
tabelaP = np.array([0, -tabelaComprimento/2, tabelaAltura])
tabela2P = np.array([bolaRaio, -tabelaComprimento/2, tabelaAltura])

######## Funções ########
def vectorModulus(V):
    return np.sqrt(np.dot(V, V))

def distancia(V1, V2):
    return np.sqrt(np.sum((V2-V1)**2))

def vectorNormalize(V):
    return V/vectorModulus(V)

def angleBetweenVectors(A, B):
    return np.arccos(np.clip(np.dot(vectorNormalize(A), vectorNormalize(B)), -1.0, 1.0))

def colisor(P, Vx, Dots, Angles, maxDotsX):
    if P[0] <= maxDotsX+bolaRaio and P[0] >= bolaRaio+1.5*deltaT*Vx:
        grid = coords2Grid(P)
        imin = max(0, int(grid[0]-np.round(bolaRaio/squareComprimento)-1))
        imax = min(gridSizeX, int(grid[0]+np.round(bolaRaio/squareComprimento)+2))
        jmin = max(0, int(grid[1]-np.round(bolaRaio/squareLargura)-1))
        jmax = min(gridSizeY, int(grid[1]+np.round(bolaRaio/squareLargura)+2))
        for i in range(imin, imax):
            for j in range(jmin, jmax):
                vectorD = P - Dots[i][j]
                if vectorModulus(vectorD) <= bolaRaio:
                    if (vectorD[0] - bolaRaio) >= 1.5*deltaT*Vx:
                        alphaLimit = np.arctan(squareLargura/(np.sqrt(2)*
                            (np.sqrt(bolaRaio**2-squareLargura**2/2) - abs(Vx)*deltaT)))
                        alpha = angleBetweenVectors(vectorD, np.array([1,0,0]))
                        #print("\n {}, {}".format(i, j))
                        #print("vectorD {}".format(vectorD))
                        #print("alphaLimit {}, alpha {}".format(np.degrees(alphaLimit), np.degrees(alpha)))
                        if alpha <= alphaLimit:
                            #print("entrou")
                            #print("vectorD - bolaRaio = {}".format(vectorD[0]-bolaRaio))
                            #print("deltaT*Vx = {}".format(deltaT*Vx))
                            return np.array([i,j])
    return False

def shoot(P0, V0, timeRange=2, rebatida=False, Dots=None, Angles=None, maxDotsX=None, grafico=False, size=1, completo=False):
    grid = 0
    X = []
    Y = []
    Z = []
    V = []
    t0 = 0
    if not rebatida and not completo:
        t0 = (maxDotsX + bolaRaio + 50 - P0[0])/V0[0]
        #print(t0)
    for t in np.arange(t0, timeRange, deltaT):
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
                    return calcErro([X[-1],Y[-1],Z[-1]], aroP, grafico=grafico)
            else:
                test = colisor((X[-1], Y[-1], Z[-1]), V0[0], Dots, Angles, maxDotsX)
                if test is not False and test[0] in range(gridSizeX) and test[1] in range(gridSizeY): #### Colidiu ####
                    #print("tempo para bater na tabela = {}".format(t))
                    grid = test
                    VF = np.array(V[-1])
                    PF = np.array([X[-1], Y[-1], Z[-1]])
                    if grafico:
                        ax.scatter3D(X, Y, Z, s=size)
                    return True, grid, PF, VF
    if grafico:
        ax.scatter3D(X, Y, Z, s=size)
    return False, False, False, False #### Não Colidiu / Não Detectou ####

def refletor(P, V, theta, phi, grafico = False):
    normal = calcNormal(theta, phi)
    if grafico:    ### Desenha a Normal do bloco colidido ###
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

def calcErro(P, Target, grafico=False):
    deltaX = Target[0] - P[0]
    deltaY = Target[1] - P[1]
    if grafico:
        ax.quiver(P[0], P[1], P[2], # <-- starting point of vector
                deltaX, deltaY, 0, # <-- directions of vector
                color = 'black', alpha = .2, lw = 2)
    return np.sqrt(deltaX**2 + deltaY**2)

def grid2Coords(i, j):
    y = i * tabelaComprimento/(gridSizeX-1)/2
    z = j * tabelaLargura/(gridSizeY-1) + tabelaAltura
    x = 0
    return np.array([x, y, z])

def coords2Grid(P):
    y = P[1]
    z = P[2] - focoAltura - tabelaAltura
    gridX = y*2*(gridSizeX-1) / tabelaComprimento
    gridY = (z/tabelaLargura + 1/2)*(gridSizeY-1)
    return np.array([round(gridX), round(gridY)])

def grid2SquareCoords(gridX, gridY):
    global squareComprimento, squareLargura
    P = grid2Coords(gridX, gridY)
    x = P[1] - squareComprimento/2
    y = P[2] - squareLargura/2
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
        fc = shoot(PF, refletor(PF, VF, theta, c, grafico = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
        fd = shoot(PF, refletor(PF, VF, theta, d, grafico = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
    else:
        fc = shoot(PF, refletor(PF, VF, c, phi, grafico = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
        fd = shoot(PF, refletor(PF, VF, d, phi, grafico = False),
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
                fc = shoot(PF, refletor(PF, VF, theta, c, grafico = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
            else:
                fc = shoot(PF, refletor(PF, VF, c, phi, grafico = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
            yc = fc
        else:
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            if var == "phi":
                fd = shoot(PF, refletor(PF, VF, theta, d, grafico = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
            else:
                fd = shoot(PF, refletor(PF, VF, d, phi, grafico = False),
                      rebatidaTimeRange, rebatida=True, grafico=False)
            yd = fd
    if yc < yd:
        #return (a, d)
        return (a + d)/2
    else:
        #return (c, b)
        return (c + b)/2

def preparaPlot(visao=1, mostraQuadra=False):
    ######### Configura o PyPlot ##########
    global fig, ax
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
    aro = plt.Circle(aroP, aroRaio, color='r', fill=0)
    mosca = plt.Circle(aroP, moscaRaio, color='blue', fill=1)
    tabela = plt.Rectangle((tabelaP[1],tabelaP[2]),
                           tabelaComprimento, tabelaLargura, color='w', fill=0)
    tabela2 = plt.Rectangle((tabela2P[1],tabela2P[2]),
                           tabelaComprimento, tabelaLargura, color='b', fill=0, ls="--")
    ax.add_patch(aro)
    ax.add_patch(mosca)
    ax.add_patch(tabela)
    ax.add_patch(tabela2)
    art3d.pathpatch_2d_to_3d(aro, z=aroP[2], zdir="z")
    art3d.pathpatch_2d_to_3d(mosca, z=aroP[2], zdir="z")
    art3d.pathpatch_2d_to_3d(tabela, z=0, zdir="x")
    art3d.pathpatch_2d_to_3d(tabela2, z=bolaRaio, zdir="x")
    if mostraQuadra:
        quadra = plt.Rectangle((aroP[0]-1575,-quadraLargura/2),
                           quadraComprimento, quadraLargura, color='black', fill=0)
        garrafao = plt.Rectangle((aroP[0]-1575,-garrafaoLargura/2),
                           garrafaoComprimento, garrafaoLargura, color='black', fill=0)
        garrafaoCirculo = plt.Circle((aroP[0]-1575+garrafaoComprimento,0),
                           garrafaoCirculoDiametro/2, color='black', fill=0)
        ax.add_patch(quadra)
        ax.add_patch(garrafao)
        ax.add_patch(garrafaoCirculo)
        art3d.pathpatch_2d_to_3d(quadra, z=0, zdir="z")
        art3d.pathpatch_2d_to_3d(garrafao, z=0, zdir="z")
        art3d.pathpatch_2d_to_3d(garrafaoCirculo, z=0, zdir="z")
    return ax

def plotNormal(i, j, Angles, x=0, scale=bolaRaio/4):
    P = grid2Coords(i, j)
    normal = scale*calcNormal(Angles[i][j][0],Angles[i][j][1])
    ax.quiver(P[0]+x, P[1], P[2], # <-- starting point of vector
    normal[0], normal[1], normal[2], # <-- directions of vector
    color = 'red', alpha = 0.5, lw = 2)
    
def plotSquare(i, j, Dots, projeta=True):
    P = grid2SquareCoords(i, j)
    gridSquare = plt.Rectangle(P, squareComprimento, squareLargura, color='blue', fill=1)
    ax.add_patch(gridSquare)
    z = Dots[i][j][0]
    if projeta:
        z += bolaRaio
    art3d.pathpatch_2d_to_3d(gridSquare, z=z, zdir="x")

def plotTabela(Dots, Angles, setas=False, scale=bolaRaio/4):
    X, Y, Z = [], [], []
    for i in range(gridSizeX):
        for j in range(gridSizeY):
            X.append(Dots[i][j][0])
            Y.append(Dots[i][j][1])
            Z.append(Dots[i][j][2])
            if setas:
                plotNormal(i, j, Angles, Dots[i][j][0], scale=scale)
    if not setas:
        ax.scatter3D(X, Y, Z, s=5)