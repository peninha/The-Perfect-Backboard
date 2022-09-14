# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 22:14:14 2022

@author: millo
"""

import numpy as np
import matplotlib.pyplot as plt

Grid2 = np.delete(Grid, np.arange(45), axis=0)

Theta = Grid2[:,:,0]
Theta2 = Theta.copy()
Theta2[:,:2] = np.nan
gxx, gyy = np.gradient(Theta2[:,2:])
Theta2[:,1] =  Theta2[:,2] - gyy[:,0]

Theta = Grid2[:,:,0]
Theta2 = Theta.copy()
Theta2[:,:1] = np.nan
gxx, gyy = np.gradient(Theta2[:,1:])
Theta2[:,0] =  Theta2[:,1] - gyy[:,0]

Phi = Grid2[:,:,1]
Phi2 = Phi.copy()
Phi2[:,:2] = np.nan
gxx, gyy = np.gradient(Phi2[:,2:])
_, ggyy = np.gradient(gyy)
Phi2[:,1] =  Phi2[:,2] - gyy[:,0] + gyy[:,0]/2

Phi = Grid2[:,:,1]
Phi2 = Phi.copy()
Phi2[:,:1] = np.nan
gxx, gyy = np.gradient(Phi2[:,1:])
_, ggyy = np.gradient(gyy)
Phi2[:,0] =  Phi2[:,1] - gyy[:,0] + gyy[:,0]/2
