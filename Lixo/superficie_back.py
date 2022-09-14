# -*- coding: utf-8 -*-
"""
Gerador de superfície para uma Tabela de Basquete que sempre acerta na cesta

Gera a nova superfície a partir dos ângulos Phi e Theta calculados

@author: Pena
"""

from funcoes import *

Angles = np.load('data/Angles2.npy')

#Angles = Angles[:,5:,:] #corta os pontos falsos

gridSizeX = len(Angles[:])
gridSizeY = len(Angles[:][0])

dotSize = 1
visao = 3
ax = preparaPlot(visao)
ax.view_init(elev=0., azim=0.01)

Theta = Angles[:,:,0]
Z = np.zeros((gridSizeX,gridSizeY))
Y = np.zeros((gridSizeX,gridSizeY))
Xx = np.zeros((gridSizeX,gridSizeY))
for j in range(gridSizeY):
    for i in range(gridSizeX):
        P = grid2Coords(i, j)
        Z[i][j] = P[2]
        Y[i][j] = P[1]
        if i==0:
            continue
        Xx[i][j] = Xx[i-1][j] + squareComprimento*(np.tan(Theta[i-1][j]) + np.tan(Theta[i][j]))/2
    
Phi = Angles[:,:,1]
Z = np.zeros((gridSizeX,gridSizeY))
Y = np.zeros((gridSizeX,gridSizeY))
Xy = np.zeros((gridSizeX,gridSizeY))
for i in range(gridSizeX):
    #for j in range(30,gridSizeY):
    for j in range(25,gridSizeY):
        #P = grid2Coords(i, j)
        P = grid2Coords(i, j)
        Z[i][j] = P[2]
        Y[i][j] = P[1]
        #if j==30 and i==0:
        if j==25 and i==0:
            continue
        Xy[i][j] = Xy[i][j-1] - squareLargura*(np.tan(Phi[i][j-1]) + np.tan(Phi[i][j]))/2
        
        k = 55-j
        #P = grid2Coords(i, k)
        P = grid2Coords(i, k)
        Z[i][k] = P[2]
        Y[i][k] = P[1]
        Xy[i][k] = Xy[i][k+1] + squareLargura*(np.tan(Phi[i][k+1]) + np.tan(Phi[i][k]))/2
X = Xy + Xx

Normals = np.zeros((gridSizeX,gridSizeY,3))
for i in range(gridSizeX):
    for j in range(gridSizeY):
        Normals[i][j] = calcNormal(Theta[i][j], Phi[i][j])
        #plotNormal(i, j, Normals, x = X[i][j])

Dots1 = np.zeros((gridSizeX,gridSizeY,3))
Dots1[:,:,0] = X
Dots1[:,:,1] = Y
Dots1[:,:,2] = Z

Angles2 = np.zeros((gridSizeX,gridSizeY,2))
gxx, gyy = np.gradient(X)
Angles2[:,:,0] = np.arctan2(gxx,10)
Angles2[:,:,1] = -np.arctan2(gyy,10)

ax.scatter3D(Xx, Y, Z, s=dotSize)
ax.scatter3D(Xy, Y, Z, s=dotSize)
#ax.scatter3D(X, Y, Z, s=dotSize)
plt.show()