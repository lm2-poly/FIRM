% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [G, loc2glob, coor_loc, area, l23, l31, l12]=local_frame(coor)

% This function returns the local, side-aligned, frame of a triangle in 3D space
%
% input:
%
% coor(3x3)                 nodal coordinates [V1, V2, V3] in global frame
%                           coor(:,i) contains the coordinates of node V_i
%
% output:
%
% G(3x1)                    centroid of the triangle, in global frame
%
% loc2glob(3x3)             local unit vectors e,h,n, in global frame
%                           first column:  e unit vector (aligned to V2-V1)
%                           second column: h unit vector (h = n x e)
%                           third column:  n unit vector (normal to triangle, oriented as (V2-V1)x(V3-V1)
%
% coor_loc(3x3)             nodal coordinates [V1,V2,V3] in local frame (G,e,h,n)
%                           the third row of coor_loc is zero
%                           coor_loc(:,i) contains the coordinates of node V_i in local frame
%
% l23,l31,l12               lengths of the triangle sides, respectively opposites to nodes V1, V2, V3
%
% area                      triangle area
%
% remark on coordinate transformation
%   loc2glob * (local components of a vector) = global components of that vector
%   loc2glob'* (global components of a vector) = local components of that vector
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

% centroid of the triangle in global frame
% equation (11)
G=sum(coor,2)/3;

% local frame
% equation (14)
% e is aligned to V12=V2-V1
V12=coor(:,2)-coor(:,1);
l12=norm(V12);
e=V12/l12;
% n is the unit vector oriented as (V2-V1)x(V3-V1)
V13=coor(:,3)-coor(:,1);
v=cross(V12,V13);
norm_v=norm(v);
n=v/norm_v;
% area is half of the norm of v
area=norm_v/2;
% h = n x e
h=cross(n,e);

% loc2glob contains columnwise the local unit vectors e,h,n in global frame
loc2glob=[e,h,n];

% nodal coordinates [V1,V2,V3] in local frame (G,e,h,n)
a=V13'*e;
b=V13'*h;
% origin at centroid
x=[0,l12,a]-(l12+a)/3;
y=[0,0,b]-b/3;
% local z is zero
z=zeros(1,3);
% assemble into a 3x3 matrix
coor_loc=[x;y;z];

% side lengths l23 and l31
bsq=b^2;
l23=sqrt((l12-a)^2+bsq);
l31=sqrt(a^2+bsq);

end
