#Meshes a parallel slits kirigami pattern on a rectangle based on its dimensions using GMSH (also called Ribbon kirigami)
#Based on Marzin et al. (2022)'s description
#Author: Danick Lamoureux
#Creation date: 2023-05-23
#See mesh_parallel_slits_kirigami.m for more information

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

def parallel_slit_rectangle_mesh(L, W, Ls, dx, dy, cut_thickness, t, mesh_size, show = False):
    Nx = int(np.round((W/(2*dx+Ls))))
    Ny = int(np.round(L/dy)-1)
    G1 = [0,0,0]
    G2 = [W,0,0]
    G3 = [W,L,0]
    G4 = [0,L,0]
    
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
    cut = cut_thickness

    wireTags = []
    #Cuts in x
    for j in range(Ny):
        if (np.mod(j,2)==0):
            for i in range(Nx+1):
                Dx = (i)*(2*dx + Ls) - 0.5*Ls
                Dy = (1+j)*dy
                k = lambda i,j: 10*(Ny*i+j) + 11
                if i == 0:
                    gmsh.model.geo.addPoint(x0+Dx+Ls/2, y0+Dy, z0, mesh_size, k(i,j))
                    gmsh.model.geo.addPoint(x0+Dx+Ls, y0+Dy, z0, mesh_size, k(i,j)+1)
                    gmsh.model.geo.addPoint(x0+Dx+Ls, y0+Dy+cut, z0, mesh_size, k(i,j)+2)
                    gmsh.model.geo.addPoint(x0+Dx+Ls/2, y0+Dy+cut, z0, mesh_size, k(i,j)+3)
                    gmsh.model.geo.addPoint(x0+Dx+Ls - cut/10, y0+Dy+cut/2, z0, mesh_size, k(i,j)+5)
                elif i == Nx:
                    gmsh.model.geo.addPoint(x0+Dx, y0+Dy, z0, mesh_size, k(i,j))
                    gmsh.model.geo.addPoint(x0+Dx+Ls/2, y0+Dy, z0, mesh_size, k(i,j)+1)
                    gmsh.model.geo.addPoint(x0+Dx+Ls/2, y0+Dy+cut, z0, mesh_size, k(i,j)+2)
                    gmsh.model.geo.addPoint(x0+Dx, y0+Dy+cut, z0, mesh_size, k(i,j)+3)
                    gmsh.model.geo.addPoint(x0+Dx+cut/10, y0+Dy+cut/2, z0, mesh_size, k(i,j)+4)
                else:
                    gmsh.model.geo.addPoint(x0+Dx, y0+Dy, z0, mesh_size, k(i,j))
                    gmsh.model.geo.addPoint(x0+Dx+Ls, y0+Dy, z0, mesh_size, k(i,j)+1)
                    gmsh.model.geo.addPoint(x0+Dx+Ls, y0+Dy+cut, z0, mesh_size, k(i,j)+2)
                    gmsh.model.geo.addPoint(x0+Dx, y0+Dy+cut, z0, mesh_size, k(i,j)+3)
                    gmsh.model.geo.addPoint(x0+Dx+Ls-cut/10, y0+Dy+cut/2, z0, mesh_size, k(i,j)+5)
                    gmsh.model.geo.addPoint(x0+Dx+cut/10, y0+Dy+cut/2, z0, mesh_size, k(i,j)+4)
                
                
                    gmsh.model.geo.addLine(k(i,j), k(i,j)+1, k(i,j)+6)
                    gmsh.model.geo.addCircleArc(k(i,j)+1, k(i,j)+5, k(i,j)+2, k(i,j)+7)
                    gmsh.model.geo.addLine(k(i,j)+2, k(i,j)+3, k(i,j)+8)
                    gmsh.model.geo.addCircleArc(k(i,j)+3, k(i,j)+4, k(i,j), k(i,j)+9)
                    gmsh.model.geo.addCurveLoop(curveTags=[k(i,j)+6, k(i,j)+7, k(i,j)+8, k(i,j)+9], tag = k(i,j)+10)
                    wireTags.append(k(i,j)+10)
        else:
            for i in range(Nx):
                Dx = (i)*(2*dx + Ls) + dx
                Dy = (1+j)*dy 
                k = lambda i,j: 10*(Ny*i+j) + 11
                gmsh.model.geo.addPoint(x0+Dx, y0+Dy, z0, mesh_size, k(i,j))
                gmsh.model.geo.addPoint(x0+Dx+Ls, y0+Dy, z0, mesh_size, k(i,j)+1)
                gmsh.model.geo.addPoint(x0+Dx+Ls, y0+Dy+cut, z0, mesh_size, k(i,j)+2)
                gmsh.model.geo.addPoint(x0+Dx, y0+Dy+cut, z0, mesh_size, k(i,j)+3)
                gmsh.model.geo.addPoint(x0+Dx+Ls-cut/10, y0+Dy+cut/2, z0, mesh_size, k(i,j)+5)
                gmsh.model.geo.addPoint(x0+Dx+cut/10, y0+Dy+cut/2, z0, mesh_size, k(i,j)+4)
            
            
                gmsh.model.geo.addLine(k(i,j), k(i,j)+1, k(i,j)+6)
                gmsh.model.geo.addCircleArc(k(i,j)+1, k(i,j)+5, k(i,j)+2, k(i,j)+7)
                gmsh.model.geo.addLine(k(i,j)+2, k(i,j)+3, k(i,j)+8)
                gmsh.model.geo.addCircleArc(k(i,j)+3, k(i,j)+4, k(i,j), k(i,j)+9)
                gmsh.model.geo.addCurveLoop(curveTags=[k(i,j)+6, k(i,j)+7, k(i,j)+8, k(i,j)+9], tag = k(i,j)+10)
                wireTags.append(k(i,j)+10)
                

    contourTags = [0, 1]
    for j in range(Ny):
        if (np.mod(j,2)==0):
            k = lambda i,j: 10*(Ny*i+j) + 11
            i = Nx
            contourTags.append(k(i,j)+1)
            contourTags.append(k(i,j))
            contourTags.append(k(i,j)+4)
            contourTags.append(k(i,j)+3)
            contourTags.append(k(i,j)+2)

    contourTags.append(2)
    contourTags.append(3)
    reversecontourTags = []
    for j in range(Ny):
        if (np.mod(j,2)==0):
            k = lambda i,j: 10*(Ny*i+j) + 11
            i=0
            reversecontourTags.append(k(i,j))
            reversecontourTags.append(k(i,j)+1)
            reversecontourTags.append(k(i,j)+5)
            reversecontourTags.append(k(i,j)+2)
            reversecontourTags.append(k(i,j)+3)

    for j in range(1,len(reversecontourTags)+1):
        contourTags.append(reversecontourTags[-j])

    startTag = gmsh.model.geo.getMaxTag(1)
    contourTags.append(0)
    contourLines = []

    i = 0 #Line counter
    k = 0 #Starting point tag counter

    newRule = False
    while k < (len(contourTags) - 1):
        if contourTags[k] == 2:
            newRule = True

        if (np.mod(i+1, 4) != 0 and not newRule) or (np.mod(i-1, 4) != 0 and newRule):
            gmsh.model.geo.addLine(contourTags[k], contourTags[k+1], startTag+i+1)
            k += 1
        else:
            gmsh.model.geo.addCircleArc(contourTags[k], contourTags[k+1], contourTags[k+2], startTag+i+1)
            k += 2
        contourLines.append(startTag+i+1)
        i += 1
        

    gmsh.model.geo.addCurveLoop(contourLines, tag = startTag+len(contourTags)+1)
    tags = [startTag+len(contourTags)+1]
    for item in wireTags:
        tags.append(item)
    gmsh.model.geo.addPlaneSurface(tags,1)
    
    dimTags = gmsh.model.getPhysicalGroups(dim = 0)
    
    gmsh.model.removePhysicalGroups(dimTags)
    
    gmsh.model.addPhysicalGroup(2, [1], 2)

    gmsh.model.geo.synchronize()

    gmsh.option.setNumber("Mesh.Smoothing", 5)

    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)
    gmsh.option.setNumber('Mesh.Points', 1)

    gmsh.option.setNumber('Mesh.MeshSizeMin', mesh_size/10)
    gmsh.option.setNumber('Mesh.MeshSizeMax', mesh_size*3)
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
        nodes[i,2] = np.random.uniform(-0.001, 0.001)*t
    
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
    Ls = float(sys.argv[3])
    dx = float(sys.argv[4])
    dy = float(sys.argv[5])
    cut_thickness = float(sys.argv[6])
    t = float(sys.argv[7])
    mesh_size = float(sys.argv[8])
    show = sys.argv[9]
    if show == "False": show = False
    else: show = True
    print("Show = "+str(show))
    nodes, simplices = parallel_slit_rectangle_mesh(L, W, Ls, dx, dy, cut_thickness, t, mesh_size, show)
    