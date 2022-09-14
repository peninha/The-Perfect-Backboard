# -*- coding: utf-8 -*-
"""
Spyder Editor

Programa por @Penadoxo
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d

def colisor(P):
    global bolaRaio
    if P[0] <= bolaRaio:
        return True
    return False

def calcV(P1, P0, deltaT):
    return (P1 - P0)/deltaT

def reflect(Q1, V1, t, deltaT):
    global P0, V0
    P0 = Q1 
    V0 = V1
    V0[0] = -V0[0]
    
def detector(z, Vz):
    global aroAltura
    if Vz <= 0:
        if z <= aroAltura:
            return True
    return False
    
def calcNormal(theta, phi):
    return np.array([np.cos(phi)*np.cos(theta), -np.cos(phi)*np.sin(theta),
                     np.sin(phi)])
    


bolaRaio = 119
aroRaio = 450/2
aroHaste = 150
aroAltura = 3000
tabelaComprimento = 900
tabelaLargura = 600
focoAltura = 300
shotAltura = 1850
shotDistancia = 4600

aroP = np.array([aroHaste + aroRaio, 0, aroAltura])
aro2P = np.array([bolaRaio - (aroP[0]-bolaRaio), 0, aroAltura])
tabelaP = np.array([0, -tabelaComprimento/2, aroAltura])
tabela2P = np.array([bolaRaio, -tabelaComprimento/2, aroAltura])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
aro = plt.Circle(aroP, aroRaio, color='r', fill=0)
aro2 = plt.Circle(aro2P, aroRaio, color='r', fill=0)
tabela = plt.Rectangle((tabelaP[1],tabelaP[2]),
                       tabelaComprimento, tabelaLargura, color='b', fill=0)
tabela2 = plt.Rectangle((tabela2P[1],tabela2P[2]),
                       tabelaComprimento, tabelaLargura, color='b', fill=0)
ax.add_patch(aro)
#ax.add_patch(aro2)
ax.add_patch(tabela)
ax.add_patch(tabela2)
art3d.pathpatch_2d_to_3d(aro, z=aroP[2], zdir="z")
#art3d.pathpatch_2d_to_3d(aro2, z=aroP[2], zdir="z")
art3d.pathpatch_2d_to_3d(tabela, z=0, zdir="x")
art3d.pathpatch_2d_to_3d(tabela2, z=bolaRaio, zdir="x")
ax.set_xlim(0, 5000)
ax.set_ylim(-2500, 2500)
ax.set_zlim(0, 5000)

#ax.set_xlim(0, 300)
#ax.set_ylim(-450, -150)
#ax.set_zlim(3000, 3300)


#Cálculo da Parábola Principal
g =  9807 #mm/s²
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
VZ0 = DeltaZ1/t1 + g*t1/2
#VZ0 = 8000
V0 = np.array([VX0, 0, VZ0])


def Px(t, P0, V0, g):
    return P0[0]+V0[0]*t

def Py(t, P0, V0, g):
    return 0

def Pz(t, P0, V0, g):
    return P0[2]+V0[2]*t-g*t**2/2

deltaT = 0.0005
timeRange = 1.6
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
         if detector(Z[-1],V[-1][2]) == True:
             #erroX
             break
    if colisor((X[-1], Y[-1], Z[-1])) and hit == False:
        reflect(np.array((X[-1], Y[-1], Z[-1])),
                V[-1], t, deltaT)
        t0 = T
        hit = True
    

ax.scatter3D(X, Y, Z, s=280)
#ax.scatter3D(0, -450, 3600, s=280)
plt.show()


Grid = [[np.array([0,0,0]) for x in range(61)] for y in range(91)]
#Grid coordenadas = Theta (horizontal), Phi (vertical), X (profundidade)

X = []
Y = []
Z = []
for i in range (91):
    y = i*10 - 450
    x1 = aroP[0]
    x2 = shotDistancia
    theta = np.arctan2(y*(x1+x2),
                       np.sqrt(y**2+x2**2)*np.sqrt(y**2+x1**2)+x2*x1-y**2)
    for j in range (5):
        Grid [i][j] = np.array([theta, 0, 0])
        normal = 30*calcNormal(theta, 0)
        x = 0
        z = j*10 + 3000
        ax.quiver(
            x, y, z, # <-- starting point of vector
            normal[0], normal[1], normal[2], # <-- directions of vector
            color = 'red', alpha = .5, lw = 1,
            )
        #X.append(x)
        #Y.append(y)
        #Z.append(z)
        
        #Grid [i][j] = np.degrees(theta)

plt.show()
    

