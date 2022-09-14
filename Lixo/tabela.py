# -*- coding: utf-8 -*-
"""
Spyder Editor

Programa por @Penadoxo
"""

import math
import numpy as np
import matplotlib.pyplot as plt



bolaRaio = 119
aroRaio = 450/2
aroHaste = 150
aroAltura = 3000
focoAltura = 300
shotAltura = 1850
shotDistancia = 5700

aroP = np.array([aroHaste + aroRaioi, 0, aroAltura])


#pontoP = np.array([400, 300])
#pontoNormal = np.array([-0.5, -1])
#pontoNormal = pontoNormal/np.linalg.norm(pontoNormal)


shotP = np.array([, -800])
V = 10



shotBounceP = np.add(pontoP, bolaRaio*pontoNormal)
shotV = np.subtract(shotBounceP, shotP)
shotV = V * shotV/np.linalg.norm(shotV)
shotBounceV = np.subtract(shotV, 2*np.dot(shotV, pontoNormal)*pontoNormal)


fig, ax = plt.subplots()
ax.set_aspect('equal', 'box')
ax.set_xlim((-1000, 1000))
ax.set_ylim((-1000, 500))
#ax.view.init(azim=120)




aro = plt.Circle((0, 0), aroRaio, color='r', fill=0)
ax.add_patch(aro)
plt.plot([pontoP[0], pontoP[0]+bolaRaio*pontoNormal[0]], [pontoP[1], pontoP[1]+bolaRaio*pontoNormal[1]], 'b-')


trajectory = np.array([shotP])
for t in range(1,100):
    trajectory = np.vstack((trajectory, shotP + shotV*t*8))

trajectory = np.vstack((trajectory, shotBounceP))
for t in range(1,100):
    trajectory = np.vstack((trajectory, shotBounceP + shotBounceV*t*8))


x, y = trajectory.T
plt.scatter(x,y, c='black', s=1400, alpha=0.1)
plt.show()
    