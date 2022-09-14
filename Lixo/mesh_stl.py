# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 19:14:46 2022

@author: Pena
"""

import numpy
from stl import mesh

# Using an existing stl file:
your_mesh = mesh.Mesh.from_file('Output/bpa_mesh2.stl')

# Or creating a new mesh (make sure not to overwrite the `mesh` import by
# naming it `mesh`):
#VERTICE_COUNT = 100
#data = numpy.zeros(VERTICE_COUNT, dtype=mesh.Mesh.dtype)
#your_mesh = mesh.Mesh(data, remove_empty_areas=False)

# The mesh normals (calculated automatically)
your_mesh.normals
Normals3 = your_mesh.normals
V0 = your_mesh.v0
V1 = your_mesh.v1
V2 = your_mesh.v2

Normals4 = Normals3/Normals3.sum(axis=1,keepdims=1)

Nquadrado = Normals3*Normals3

Normals4 = Normals3/np.sqrt(Nquadrado.sum(axis=1,keepdims=1))

# The mesh vectors
#your_mesh.v0, your_mesh.v1, your_mesh.v2
# Accessing individual points (concatenation of v0, v1 and v2 in triplets)
#assert (your_mesh.points[0][0:3] == your_mesh.v0[0]).all()
#assert (your_mesh.points[0][3:6] == your_mesh.v1[0]).all()
#assert (your_mesh.points[0][6:9] == your_mesh.v2[0]).all()
#assert (your_mesh.points[1][0:3] == your_mesh.v0[1]).all()

#your_mesh.save('new_stl_file.stl')