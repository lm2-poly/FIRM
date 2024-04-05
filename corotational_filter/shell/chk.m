% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [ ... % tensors needed to filter out the rigid rotation chk_R
    chk_R, chk_G, chk_P, ...
    ... % core nodal displacements and rotation tensors
    bar_u_i, bar_R_i]=chk( ...
    ... % nodal coordinates of reference triangle
    coor, ...
    ... % nodal displacements and rotation tensors after filtering out the rigid rototranslation t,hat_R
    hat_u_i, hat_R_i, ...
    ... % area and side lengths of reference triangle
    area, l23, l31, l12, ...
    ... % area and side lengths of deformed triangle
    area_def, l23_def, l31_def, l12_def)

% This functions computes and filters out the rigid rotation chk_R
%
% computations are performed using coordinates in the element reference frame (e,h,n)
% (i.e., e is parallel to 1-2 side of reference triangle,
% n is normal to reference triangle
% the origin is at the centroid
%
% input:
%
% coor(3x3)                 nodal coordinates of reference triangle (in element reference frame)
%                           the third row (z-coordinate) is zero
%                           coor(:,i) contains the coordinates of node V_i
%
% hat_u_i(3x3), hat_R_i(3x3x3) nodal displacements and rotation tensors obtained
%                           after filtering out the rigid rototranslation t,hat_R
%                           hat_u_i(:,i) and hat_R_i(:,:,3) are relevant to node i
%
% area, l23, l31, l12       area and side lengths of reference triangle
%
% area_def, l23_def, l31_def, l12_def   area and side lengths of deformed triangle
%
% output:
%
% chk_R, chk_G, chk_P       tensors needed to filter out the rigid rotation chk_R
%   chk_R(3x3)              tensor chk_R: equation (19)
%                           delta_chk_R defines delta_chk_theta=ax(chk_Rt*delta_chk_R)
%   chk_G(3x18)             delta_chk_theta = chk_G*shp_chk_Rt*delta_hat_a
%   chk_P(18x18)            delta_chk_a = chk_P*shp_chk_Rt*delta_hat_a
%
%
% bar_u_i(3x3), bar_R_i(3x3x3) core nodal displacements and rotation tensors,
%                           obtained after filtering out also the rigid rotation chk_R
%                           bar_u_i(:,i) and bar_R_i(:,:,3) are relevant to node i
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

% number of element nodes
n_el_n = 3;

% number of element degrees of freedom
n_el_dof = 6*n_el_n;

% indices of translational and rotational dofs
[ind_transl, ind_rot]=indices(n_el_n);

% nodal coordinates of reference triangle in the element reference frame
x=coor(1,:);
y=coor(2,:);

% computation of chk_R
% chk_R is the rotation tensor obtained from the polar decomposition
% of the gradient of the in-plane homogeneous deformation
% mapping the reference triangle onto the deformed triangle rotated back by hat_Rt

% intrinsic computation
% it is simplified by the observation that the reference and
% deformed-rotated-back triangles have their 1-2 sides aligned
% equation (20)
% c1,c2,c3 scalars
c1=(l31^2+l12^2-l23^2)*l23_def^2;
c2=(l12^2+l23^2-l31^2)*l31_def^2;
c3=(l23^2+l31^2-l12^2)*l12_def^2;
% trace of the deformation tensor Uh
tr_Uh=sqrt( 2*area_def/area + (c1+c2+c3)/(8*area^2) );
% angle alpha, characterizing the rotation tensor chk_R
d=area*l12*l12_def*tr_Uh;
cos_alpha=(area*l12_def^2+area_def*l12^2)/d;
sin_alpha=(l12^2*(l23_def^2-l31_def^2)+l12_def^2*(l31^2-l23^2))/(4*d);

% rotation tensor chk_R
% equation (19)
chk_R=[[cos_alpha;sin_alpha;0], [-sin_alpha;cos_alpha;0], [0;0;1]];

% displacement bar_u_i, obtained after filtering out the rigid rotation chk_R
% equation (24), the origin is at the centroid
bar_u_i=chk_R'*(coor+hat_u_i)-coor;

% rotation tensor bar_R_i at each node, obtained after filtering out the rigid rotation chk_R
% equation (26)
bar_R_i=zeros(3,3,n_el_n);
for i=1:n_el_n
    bar_R_i(:,:,i)=chk_R'*hat_R_i(:,:,i);
end

% tensors chk_Y(3xn_el_dof) and chk_Csi(3x3), defining the tensor chk_G (3xn_el_dof)
% chk_Y is independent from nodal parameters
% chk_Y
% equation (62)
chk_Y=zeros(3,n_el_dof);
chk_Y(3,:) = [ ...
    (x(2)-x(3)), (y(2)-y(3)), 0,   ...
    0,           0,           0,   ...
    (x(3)-x(1)), (y(3)-y(1)), 0,   ...
    0,           0,           0,   ...
    (x(1)-x(2)), (y(1)-y(2)), 0,   ...
    0,           0,           0,   ...
    ]/(2*area);

% chk_Csi
% equation (61)
chk_Csi=eye(3)/tr_Uh;

% tensor chk_G (3xn_el_dof), defined by:
% delta_chk_theta = chk_G*shp_chk_Rt*delta_hat_a, equation (55)
% equation (60)
chk_G=chk_Csi*chk_Y;

% tensor chk_A (n_el_dofx3), stack(-chk_A_i,I),
% with chk_A_i=spin(V_i+chk_u_i-G)
% stack of spins of position vectors
% equation (69)
chk_A=zeros(n_el_dof,3);
% nodal coordinates of bar triangle (i.e., triangle deformed by Uh)
% in element reference frame
coor_bar=coor+bar_u_i;
for i=1:n_el_n
    ind_transl_i=ind_transl(:,i);
    ind_rot_i=ind_rot(:,i);
    chk_A(ind_transl_i,:)=-spin(coor_bar(:,i));
    chk_A(ind_rot_i,:)=eye(3);
end

% projector chk_P (n_el_dofxn_el_dof), defined by:
% delta_chk_a = chk_P*shp_chk_Rt*delta_hat_a, equation (67)
% equation (68)
chk_P=eye(n_el_dof)-chk_A*chk_G;

end
