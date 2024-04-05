#Meshes a perforated plate based on its dimensions using GMSH. Based on Jin et al. (2020)'s description
#Y. Jin, J.-T. Kim, S. Cheng, O. Barry, and L. P. Chamorro, On the distinct drag, reconfiguration and wake of perforated
#structures, Journal of Fluid Mechanics 10.1017/jfm.2020.98 (2020)
#Author: Danick Lamoureux
#Creation date: 2023-05-23
#See mesh_perforated_plate.m for more information

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

def perforated_plate(L, W, a, w, mesh_size, show = False):
    Nx = int(np.floor(((L)/(a+w))))
    Ny = int(np.floor((W)/(a+w)))
    G1 = [0,0,0]
    G2 = [L,0,0]
    G3 = [L,W,0]
    G4 = [0,W,0]
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("model")
    
    geom = gmsh.model.geo
    
    gmsh.model.geo.addPoint(G1[0], G1[1], G1[2], mesh_size, 0)
    gmsh.model.geo.addPoint(G2[0], G2[1], G2[2], mesh_size, 1)
    gmsh.model.geo.addPoint(G3[0], G3[1], G3[2], mesh_size, 2)
    gmsh.model.geo.addPoint(G4[0], G4[1], G4[2], mesh_size, 3)
    
    x0 = G1[0]
    y0 = G1[1]
    z0 = G1[2]
    cut = a

    wireTags = []
    #Cuts in x
    for j in range(Ny):
        for i in range(Nx):
            Dx = (i)*(a+w) + w/2
            Dy = (j)*(a+w) + w/2
            k = lambda i,j: 8*(Ny*i+j) + 9
            gmsh.model.geo.addPoint(x0+Dx, y0+Dy, z0, mesh_size, k(i,j))
            gmsh.model.geo.addPoint(x0+Dx+cut, y0+Dy, z0, mesh_size, k(i,j)+1)
            gmsh.model.geo.addPoint(x0+Dx+cut, y0+Dy+cut, z0, mesh_size, k(i,j)+2)
            gmsh.model.geo.addPoint(x0+Dx, y0+Dy+cut, z0, mesh_size, k(i,j)+3)
            for m in range(3):
                gmsh.model.geo.addLine(k(i,j)+m, k(i,j)+m+1, k(i,j)+4+m)
            gmsh.model.geo.addLine(k(i,j)+3, k(i,j), k(i,j)+7)
            gmsh.model.geo.addCurveLoop(curveTags=[k(i,j)+4, k(i,j)+5, k(i,j)+6, k(i,j)+7], tag = k(i,j)+8)
            wireTags.append(k(i,j)+8)
    contourTags = [0, 1, 2, 3]

    startTag = gmsh.model.geo.getMaxTag(1)
    contourTags.append(0)
    contourLines = []
    for i in range(len(contourTags)-1):
        gmsh.model.geo.addLine(contourTags[i], contourTags[i+1], startTag+i+1)
        contourLines.append(startTag+i+1)

    gmsh.model.geo.addCurveLoop(contourLines, tag = startTag+len(contourTags)+1)
    tags = [startTag+len(contourTags)+1]
    for item in wireTags:
        tags.append(item)
    gmsh.model.geo.addPlaneSurface(tags,1)
    
    dimTags = gmsh.model.getPhysicalGroups(dim = 0)
    
    gmsh.model.removePhysicalGroups(dimTags)
    
    gmsh.model.addPhysicalGroup(2, [1], 2)

    gmsh.model.geo.synchronize()
    # Plotting the mesh
    if show:
        gmsh.fltk.run()
    gmsh.option.setNumber("Mesh.Smoothing", 5)

    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
    gmsh.option.setNumber('Mesh.Points', 1)

    gmsh.option.setNumber('Mesh.MeshSizeMin', mesh_size)
    gmsh.option.setNumber('Mesh.MeshSizeMax', mesh_size)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.Format", 31)
    
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
    L = float(sys.argv[1])
    W = float(sys.argv[2])
    a = float(sys.argv[3])
    w = float(sys.argv[4])
    mesh_size = float(sys.argv[5])
    show = sys.argv[6]
    if show == "False": show = False
    else: show = True
    print("Show = "+str(show))
    nodes, simplices = perforated_plate(L, W, a, w, mesh_size, show = True)
    