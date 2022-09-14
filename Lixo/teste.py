# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 17:02:53 2022

@author: Pena
"""

from funcoes import *

import numpy as np
import open3d as o3d

input_path="Sample/"
output_path="Output/"
dataname="sample.xyz"
#point_cloud= np.loadtxt(input_path+dataname,skiprows=1)


Dots = np.load('data/Dots1.npy')
Angles = np.load('data/Angles1.npy')
Angles2 = np.load('data/Angles1_CalcPorGrad.npy')
Normals = np.load('data/Normals1.npy')

Normals2 = np.zeros((gridSizeX,gridSizeY,3))
for i in range(gridSizeX):
    for j in range(gridSizeY):
        Normals2[i][j] = calcNormal(Angles2[i][j][0], Angles2[i][j][1])
        #plotNormal(i, j, Normals, x = X[i][j])


Dots = Dots.reshape([2806,3])
Angles = Angles.reshape([2806,2])
Angles2 = Angles2.reshape([2806,2])
Normals = Normals.reshape([2806,3])
Normals2 = Normals2.reshape([2806,3])

pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(Dots)
#pcd.colors = o3d.utility.Vector3dVector(point_cloud[:,3:6]/255)
pcd.normals = o3d.utility.Vector3dVector(Normals)

distances = pcd.compute_nearest_neighbor_distance()
avg_dist = np.mean(distances)
radius = 3 * avg_dist
bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(
    pcd,o3d.utility.DoubleVector([radius, radius * 2]))
#dec_mesh = bpa_mesh.simplify_quadric_decimation(100000)

#dec_mesh.remove_degenerate_triangles()
#dec_mesh.remove_duplicated_triangles()
#dec_mesh.remove_duplicated_vertices()
#dec_mesh.remove_non_manifold_edges()

bpa_mesh.remove_degenerate_triangles()
bpa_mesh.remove_duplicated_triangles()
bpa_mesh.remove_duplicated_vertices()
bpa_mesh.remove_non_manifold_edges()

"""
poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
    pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
bbox = pcd.get_axis_aligned_bounding_box()
p_mesh_crop = poisson_mesh.crop(bbox)
"""

o3d.io.write_triangle_mesh(output_path+"bpa_mesh2.stl", bpa_mesh)
#o3d.io.write_triangle_mesh(output_path+"p_mesh_c.stl", p_mesh_crop)

#o3d.visualization.draw_geometries([pcd])
o3d.visualization.draw_geometries([bpa_mesh])
#o3d.visualization.draw_geometries([dec_mesh])
#o3d.visualization.draw_geometries([p_mesh_crop])

