#Meshes a disk with a single cut based on its dimension using GMSH
#Based on Schouveiler & Boudaoud (2006)'s description
#L. Schouveiler and A. Boudaoud, The rolling up of sheets in a steady flow, J. Fluid Mech. 10.1017/S0022112006000851
#(2006).
#Author: Danick Lamoureux
#Creation date: 2023-05-20
#See mesh_singleradius_disk.m for more information

#Copyright (C) 2024 Danick Lamoureux

#This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or any later version.

#This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with 
# this program. If not, see https://www.gnu.org/licenses/.

import gmsh
import numpy as np
import os

def mesh_cut_disk(R, tcut, mesh_size, show = False):
    G1 = [0, 0, 0]
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("model")
    
    geom = gmsh.model.geo
    
    p = lambda theta: [R*np.cos(theta), R*np.sin(theta)]
    G2 = [p(0)[0], p(0)[1], 0.0]
    G3 = [p(2*np.pi/3)[0], p(2*np.pi/3)[1], 0.0]
    G4 = [p(4*np.pi/3)[0], p(4*np.pi/3)[1], 0.0]
    G5 = [p(np.arcsin(-tcut/R))[0], p(np.arcsin(-tcut/R))[1], -tcut]
    G6 = [0, -tcut, 0]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("model")

    geom = gmsh.model.geo

    gmsh.model.geo.addPoint(0.0, 0.0, 0.0, mesh_size/2, 0)
    gmsh.model.geo.addPoint(G2[0], G2[1], G2[2], mesh_size, 1)
    gmsh.model.geo.addPoint(G3[0], G3[1], G3[2], mesh_size, 2)
    gmsh.model.geo.addPoint(G4[0], G4[1], G4[2], mesh_size, 3)
    gmsh.model.geo.addPoint(G5[0], G5[1], G5[2], mesh_size, 4)
    gmsh.model.geo.addPoint(G6[0], G6[1], G6[2], mesh_size, 5)


    gmsh.model.geo.addCircleArc(1,0,2,5)
    gmsh.model.geo.addCircleArc(2,0,3,6)
    gmsh.model.geo.addCircleArc(3,0,4,7)

    gmsh.model.geo.addLine(0, 1, 8)
    gmsh.model.geo.addLine(0, 2, 9)
    gmsh.model.geo.addLine(0, 3, 10)
    gmsh.model.geo.addLine(5, 4, 11)
    gmsh.model.geo.addLine(0, 5, 12)

    gmsh.model.geo.addCurveLoop(curveTags=[5], tag = 5)
    gmsh.model.geo.addCurveLoop(curveTags=[6], tag = 6)
    gmsh.model.geo.addCurveLoop(curveTags=[7], tag = 7)
    gmsh.model.geo.addCurveLoop(curveTags=[8], tag = 8)
    gmsh.model.geo.addCurveLoop(curveTags=[9], tag = 9)
    gmsh.model.geo.addCurveLoop(curveTags=[10], tag = 10)
    gmsh.model.geo.addCurveLoop(curveTags=[11], tag = 11)
    gmsh.model.geo.addCurveLoop(curveTags=[12], tag = 12)
    gmsh.model.geo.addPlaneSurface([8, 5, -9],12)   
    gmsh.model.geo.addPlaneSurface([9, 6, -10],13)   
    gmsh.model.geo.addPlaneSurface([10, 7, -11, -12],14)    
    
    gmsh.model.addPhysicalGroup(2, [0], 2)
    
    gmsh.model.geo.synchronize()
    
    gmsh.option.setNumber("Mesh.Smoothing", 5)

    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
    gmsh.option.setNumber('Mesh.Points', 1)

    gmsh.option.setNumber('Mesh.MeshSizeMin', mesh_size/10)
    gmsh.option.setNumber('Mesh.MeshSizeMax', mesh_size)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.Format", 31)
    gmsh.model.geo.synchronize()
    # Plotting the mesh
    if show:
        gmsh.fltk.run()
    
    #Getting all the nodes
    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes(includeBoundary = True)
    nNodes = int(max(nodeTags))
    nodes_coords = np.zeros([nNodes,2])
    for i in range(len(nodeTags)):
        index = int(nodeTags[i]-1)

        nodes_coords[index, 0] = coord[int(3*i)]
        nodes_coords[index, 1] = coord[int(3*i+1)]
    
    nodes = np.zeros([len(nodes_coords), 3])
    for i in range(len(nodes_coords)):
        nodes[i,0] = nodes_coords[i,0]
        nodes[i,1] = nodes_coords[i,1]
        nodes[i,2] = 0.0
    
    filename = 'Polygon.bdf'
    gmsh.write(filename)
    
    #Finding the simplices from the written file
    f = open(filename,'r')
    rows = f.readlines()
    simplices = []
    for i in range(len(rows)):
        row = rows[i].split(' ')
        row = [item for item in row if item != '']
        if row[0] == 'CTRIA3':
            new_simplice = [int(row[3]),int(row[4]),int(row[5])]
            simplices.append(new_simplice)
    simplices = np.array(simplices)
    f.close()
    
    gmsh.finalize()
    
    os.remove(filename)
    
    np.savetxt("nodes.csv", nodes)
    np.savetxt("connec.csv", simplices)
    
    return nodes, simplices

if __name__ == '__main__':
    import sys
    R = float(sys.argv[1])
    tcut = float(sys.argv[2])
    mesh_size = float(sys.argv[3])
    show = sys.argv[4]
    if show == "False": show = False
    else: show = True
    nodes, simplices = mesh_cut_disk(R, tcut, mesh_size, show)
    
    