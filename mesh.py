# -*- coding: utf-8 -*-
"""
Gerador de superfície para uma Tabela de Basquete que sempre acerta na cesta

Cria um mesh (superfície) a partir da matriz de pontos

@author: Pena
"""

from funcoes import *
import open3d as o3d

######## Configurações ##########
Dots = np.load('data/Dots1.npy')
Angles = np.load('data/Angles1.npy')
#Angles2 = np.load('data/Angles2_CalcPorGrad.npy')
Normals = np.load('data/Normals1.npy')

output_path="Output2/"

######## Adequa as Matrizes para uma única linha ########

Dots = Dots.reshape([2576,3])
Angles = Angles.reshape([2576,2])
Normals = Normals.reshape([2576,3])
#Angles2 = Angles2.reshape([2576,2])
#Normals2 = Normals2.reshape([2576,3])

######## Gera a Nuvem de Pontos no Open3D #########
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(Dots)
#pcd.colors = o3d.utility.Vector3dVector(point_cloud[:,3:6]/255)
pcd.normals = o3d.utility.Vector3dVector(Normals)

######## Método Ball-Pivoting Algorithm (BPA) #######
distances = pcd.compute_nearest_neighbor_distance()
avg_dist = np.mean(distances)
radius = 3 * avg_dist
bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(
    pcd,o3d.utility.DoubleVector([radius, radius * 2]))
#dec_mesh = bpa_mesh.simplify_quadric_decimation(100000)
bpa_mesh.remove_degenerate_triangles()
bpa_mesh.remove_duplicated_triangles()
bpa_mesh.remove_duplicated_vertices()
bpa_mesh.remove_non_manifold_edges()

o3d.io.write_triangle_mesh(output_path+"Dots0_bpa_mesh.stl", bpa_mesh)

#o3d.visualization.draw_geometries([pcd])
o3d.visualization.draw_geometries([bpa_mesh])

