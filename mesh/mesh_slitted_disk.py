#Meshes a disk with multiple cuts along its radii based on its dimension using GMSH
#Based on Gosselin et al. (2010)'s description
#F. Gosselin, E. de Langre, and B. A. Machado-Almeida, Drag reduction of flexible plates by reconfiguration, Journal of
#Fluid Mechanics 10.1017/S0022112009993673 (2010)
#Author: Danick Lamoureux
#Creation date: 2023-05-20
#See mesh_slitted_disk.m for more information

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

def mesh_slitted_disk(Re, Ri, N, tcut, mesh_size, show = False):
    G1 = [0, 0, 0]
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("model")
    
    geom = gmsh.model.geo
    
    p = lambda R, theta: [R*np.cos(theta), R*np.sin(theta)]
    gmsh.model.geo.addPoint(G1[0], G1[1], G1[2], mesh_size, 0)
    
    curves = []
    for i in range(N):
        theta1 = i*(2*np.pi/N)
        theta2 = theta1 + tcut/Ri
        G0 = p(1.01*Ri, 0.5*(theta1+theta2))
        gmsh.model.geo.addPoint(G0[0], G0[1], 0.0, 3*i+1)
        G1 = p(Re, theta1)
        gmsh.model.geo.addPoint(G1[0], G1[1], 0.0, 3*i+2)
        gmsh.model.geo.addLine(3*i+1, 3*i+2, 3*i+1)

        G2 = p(Re, theta2)
        gmsh.model.geo.addPoint(G2[0], G2[1], 0.0, 3*i+3)
        gmsh.model.geo.addLine(3*i+1, 3*i+3, 3*i+2)
        if i != 0:
            gmsh.model.geo.addCircleArc(3*i+2, 0, 3*i, 3*i+3)
        if i == N-1:
            gmsh.model.geo.addCircleArc(2, 0, 3*i+3, 3)
        lastTag = 3*i+3
    
    
    inner_circle = []
    for i in range(N):
        curves.append(-(3*i+3))
        curves.append(-(3*i+1))
        curves.append(3*i+2)
        
        theta1 = i*(2*np.pi/N)
        theta2 = theta1 + tcut/Ri
        theta3 = (i+1)*(2*np.pi/N)
        theta4 = theta3 + tcut/Ri
        
        if i != N-1:
            if i == 0:
                initialtag = lastTag+1
                prevTag = initialtag
                G0 = p(Ri, 0.5*(theta1+theta2))
                gmsh.model.geo.addPoint(G0[0], G0[1], 0.0, prevTag)
            
            G1 = p(Ri, 0.5*(theta3+theta4))
            
            gmsh.model.geo.addPoint(G1[0], G1[1], 0.0, prevTag+1)
            gmsh.model.geo.addCircleArc(prevTag, 0, prevTag+1, lastTag+1)
            prevTag += 1
        if i == N-1:
            G0 = p(Ri, 0.5*(theta1+theta2))
            gmsh.model.geo.addPoint(G1[0], G1[1], 0.0, prevTag+1)
            gmsh.model.geo.addCircleArc(initialtag, 0, prevTag+1, lastTag+1)
        inner_circle.append(lastTag+1)
        lastTag += 1
    gmsh.model.geo.addCurveLoop(curveTags=curves, tag = lastTag+1, reorient=True)
    gmsh.model.geo.addCircleArc(inner_circle[-1], 0, inner_circle[0], lastTag+3)
    inner_circle = inner_circle[0:-1]
    inner_circle.append(lastTag+3)
    gmsh.model.geo.addCurveLoop(curveTags=inner_circle, tag = lastTag+2, reorient=True)
    #gmsh.model.geo.addLine(inner_circle[-1], inner_circle[0], lastTag+3)
    #gmsh.model.geo.addCurveLoop(curveTags=inner_circle[0:-1], tag = lastTag+2, reorient=True)
    gmsh.model.geo.addPlaneSurface([lastTag+1, lastTag+2],0)
    
    gmsh.model.addPhysicalGroup(2, [0], 2)
    
    gmsh.model.geo.synchronize()
    
    gmsh.option.setNumber("Mesh.Smoothing", 5)

    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
    gmsh.option.setNumber('Mesh.Points', 1)

    gmsh.option.setNumber('Mesh.MeshSizeMin', mesh_size)
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
    Re = float(sys.argv[1])
    Ri = float(sys.argv[2])
    N = int(sys.argv[3])
    tcut = float(sys.argv[4])
    mesh_size = float(sys.argv[5])
    show = sys.argv[6]
    
    if show == "False": show = False
    else: show = True
    nodes, simplices = mesh_slitted_disk(Re, Ri, N, tcut, mesh_size, show)
    
    