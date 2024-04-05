% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [ ...
    ... % reference triangle: centroid G and element reference frame (e,h,n) (in global frame)
    G, ref2glob, ...
    ... % reference triangle: nodal coordinates (in element reference frame), area, side lengths
    coor, area, l23, l31, l12, ...
    ... % deformed triangle: centroid G_def and attatched frame (e_def,h_def,n_def) (in global frame)
    G_def, def2glob, ...
    ... % deformed triangle: nodal coordinates (in frame attatched to deformed triangle), area, side lengths
    coor_def, area_def, l23_def, l31_def, l12_def, ...
    ... % deformed triangle: attatched frame (in element reference frame)
    e_def, h_def, n_def, ...
    ... % nodal parameters (in element reference frame)
    a, ...
    ... % dead and follower loads (in element reference frame)
    d_load, f_load]=geometry( ...
    ... % nodal coordinates and nodal parameters (in global frame)
    coor_el, a_el, ...
    ... % dead and follower loads (the former in global frame; the latter in element reference frame)
    d_load_el, f_load_el)

% This function computes geometrical features of reference and deformed triangle
% (more precisely, the latter is the triangle joining corrent nodal positions);
% moreover, it prepares computations in element reference frame
%
% input:
%
% coor_el(3x3)              nodal coordinates of reference triangle (in global frame)
%                           coor_el(:,i) contains the coordinates of node V_i
%
% a_el(18x1)                nodal parameters (in global frame)
%                           a_el=(u_1;theta_1;u_2;theta_2;u_3;theta_3)
%
% d_load_el(6x3), f_load_el(6x3) dead and follower loads per unit area, linearly varying over the triangle
%                           d_load(:,i) and f_load(:,i) are load intensities at node V_i
%                           d_load is expressed in global frame; f_load is expressed in element reference frame
%
% output:
%
% G(3x1)                    centroid of reference triangle, in global frame
%
% ref2glob(3x3)             unit vectors (e,h,n) of element reference frame, in global frame
%                           first column:  e unit vector (aligned to V2-V1)
%                           second column: h unit vector (h = n x e)
%                           third column:  n unit vector (normal to triangle, oriented as (V2-V1)x(V3-V1)
%
% coor(3x3)                 nodal coordinates [V1,V2,V3] of reference triangle in element reference frame (G,e,h,n)
%                           the third row of coor is zero
%
% area                      area of reference triangle
%
% l23,l31,l12               lengths of reference triangle sides, respectively opposites to nodes V1, V2, V3
%
% G_def(3x1)                centroid of deformed triangle (joining nodes in current position), in global frame
%
% def2glob(3x3)             unit vectors (e_def,h_def,n_def) of frame attatched to deformed triangle, in global frame
%                           first column:  e_def unit vector (aligned to V2_def-V1_def)
%                           second column: h_def unit vector (h_def = n_def x e_def)
%                           third column:  n_def unit vector (normal to triangle, oriented as (V2_def-V1_def)x(V3_def-V_def1)
%
% coor_def(3x3)             nodal coordinates [V1_def,V2_def,V3_def] of deformed triangle in frame (G_def,e_def,h_def,n_def)
%
% area_def                  area of deformed triangle
%
% l23_def,l31_def,l12_def   lengths of deformed triangle sides, respectively opposites to nodes V1_def, V2_def, V3_def
%
% e_def(3x1), h_def(3x1), n_def(3x1) frame attached to deformed triangle
%
% a(18x1)                   nodal parameters (in element reference frame)
%                           a=(u_1;theta_1;u_2;theta_2;u_3;theta_3)
%
% d_load(6x3), f_load(6x3)  dead and follower loads (in element reference frame)
%
% remark on coordinate transformation
%
%   ref2glob * (components of a vector in element reference frame) = global components of that vector
%   ref2glob'* (global components of a vector) = components of that vector in element reference frame
%
%   def2glob * (components of a vector in frame attatched to deformed triangle) = global components of that vector
%   def2glob'* (global components of a vector) = components of that vector in frame attatched to deformed triangle
%
% This code is part of a Matlab toolkit distributed as supplementary material of the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528
%
% Authors' e-mail addresses:
% caselli@ing.uniroma2.it (Federica Caselli)
% bisegna@uniroma2.it (Paolo Bisegna)
%
% (C) 2010-2013 Paolo Bisegna and Federica Caselli. License: GNU General Public License (GPLv3)

% indices of translational and rotational dofs
[ind_transl, ind_rot]=indices(3);

% nodal displacements
% a_el(ind_transl) yields [[u1; v1; w1], [u2; v2; w2], [u3; v3; w3]],
% i.e., the displacement components (u,v,w) of each node 1,2,3
u_el=a_el(ind_transl);

% nodal rotation vectors
% a_el(ind_rot) yields [[q1; r1; s1], [q2; r2; s2], [q3; r3; s3]],
% i.e., the rotation vector (q,r,s) of each node 1,2,3
theta_el=a_el(ind_rot);

% element reference frame (e,h,n)
% i.e., local frame of the reference triangle
[G, ref2glob, coor, area, l23, l31, l12]=local_frame(coor_el);

% nodal coordinates of the deformed triangle
coor_def=coor_el+u_el;

% frame attatched to deformed triangle (e_def,h_def,n_def)
% i.e., local frame of the deformed triangle
[G_def, def2glob, coor_def, area_def, l23_def, l31_def, l12_def]=local_frame(coor_def);

% (e_def,h_def,n_def) expressed in element reference frame (e,h,n)
def2ref=ref2glob'*def2glob;
e_def=def2ref(:,1);
h_def=def2ref(:,2);
n_def=def2ref(:,3);

% convert to element reference frame (e,f,g)

% nodal displacements
u=ref2glob'*u_el;
% nodal rotations
theta=ref2glob'*theta_el;
% nodal parameters
a=[u; theta];
% arrange into a 18x1 vector
a=a(:);

% dead loads: they are given in global frame and now converted to element reference frame
d_load=blkdiag(ref2glob',ref2glob')*d_load_el;
% follower loads: they are already given in element reference frame
f_load=f_load_el;

end
