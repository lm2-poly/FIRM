% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [ ... % tensors needed to filter out the rigid rototranslation t,hat_R
    t, hat_R, T, hat_G, hat_P, ...
    ... % nodal displacements and rotation tensors after filtering out the rigid rototranslation t,hat_R
    hat_u_i, hat_R_i]=hat( ...
    ... % nodal coordinates of reference triangle
    coor, ...
    ... % nodal displacements and rotation tensors
    u_i, R_i, ...
    ... % area and side lengths of deformed triangle
    area_def, l23_def, l31_def, l12_def, ...
    ... % frame attached to deformed triangle
    e_def, h_def, n_def)

% This function computes and filters out the rigid rototranslation t,hat_R
%
% computations are performed in the element reference frame (e,h,n)
% (i.e., e is parallel to 1-2 side of reference triangle,
% n is normal to reference triangle)
% the origin is at the centroid
%
% input:
%
% coor(3x3)                 nodal coordinates of reference triangle in element reference frame
%                           the third row (z-coordinate) is zero
%                           coor(:,i) contains the coordinates of node V_i
%
% u_i(3x3), R_i(3x3x3)      nodal displacements and rotation tensors
%                           u_i(:,i) and R_i(:,:,3) are relevant to node i
%
% area_def, l23_def, l31_def, l12_def   area and side lengths of deformed triangle
%
% e_def(3x1), h_def(3x1), n_def(3x1) frame attached to deformed triangle
%                           (i.e., e_def parallel to 1-2 side of deformed triangle,
%                           n_def normal to deformed triangle)
%
% output:
%
% t, hat_R, T, hat_G, hat_P vector and tensors needed to filter out the rigid rototranslation t,hat_R
%   t(3x1)                  translation from centroid of reference triangle to centroid of deformed triangle
%   hat_R(3x3)              rotation tensor hat_R, mapping (e,h,n) onto (e_def,h_def,n_def)
%                           delta_hat_R defines delta_hat_theta=ax(hat_Rt*delta_hat_R)
%   T(3x18)                 delta_tau=T*shp_hat_Rt*delta_tilde_a, with delta_tau=hat_Rt*delta_t
%   hat_G(3x18)             delta_hat_theta = hat_G*shp(hat_Rt)*delta_tilde_a
%   hat_P(18x18)            delta_hat_a = hat_P*shp(hat_Rt)*delta_tilde_a
%
%
% hat_u_i(3x3), hat_R_i(3x3x3) nodal displacements and rotation tensors obtained
%                           after filtering out the rigid rototranslation t,hat_R
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

% translation from centroid of reference triangle to centroid of deformed triangle
% equation (11)
t=sum(u_i,2)/n_el_n;

% rotation tensor hat_R, mapping (e,h,n) onto (e_def,h_def,n_def)
% delta_hat_R defines delta_hat_theta=ax(hat_Rt*delta_hat_R)
% hat_R=e_def*e'+h_def*h'+n_def*n';
% equation (17)
hat_R=[e_def, h_def, n_def];

% in global coordinates this computation would be:
% ref2glob=[e,h,n];
% def2glob=[e_def,h_def,n_def];
% hat_R=ref2glob'*def2glob;
% t=ref2glob'*(G_def-G)

% displacement hat_u_i,
% obtained after the rigid motion (t,hat_R) has been filtered out
% of course, the third row of hat_u_i is zero
% equation (23), the origin is at the centroid
hat_u_i=hat_R'*(coor+u_i-repmat(t,1,n_el_n))-coor;

% rotation tensor hat_R_i at each node,
% obtained after filtering out the rigid rototranslation t,hat_R
% equation (25)
hat_R_i=zeros(3,3,n_el_n);
for i=1:n_el_n
    hat_R_i(:,:,i)=hat_R'*R_i(:,:,i);
end

% tensor T (3xn_el_dof), such that:
% delta_tau=T*shp_hat_Rt*delta_tilde_a (equation (40)),
% with delta_tau=hat_Rt*delta_t (equation (38))
% equation (42)
T=zeros(3,n_el_dof);
for i=1:n_el_n
    ind_transl_i=ind_transl(:,i);
    T(:,ind_transl_i)=eye(3)/n_el_n;
end

% tensors hat_Y(3xn_el_dof) and hat_Csi(3x3), defining the tensor hat_G(3xn_el_dof)
% hat_Y is independent from nodal parameters

% hat_Csi
% equation (44)
hat_Csi=[-l12_def/(4*area_def) , (l31_def^2-l23_def^2)/(4*l12_def*area_def) , 0; ...
    0 , 1/l12_def , 0; ...
    0 , 0 , 1/l12_def];

% hat_Y
% equation (45)
hat_Y=[ ...
    [ ...
    0 ,  0 ,  1; ...
    0 ,  0 ,  1; ...
    0 , -1 ,  0 ], ...
    zeros(3), ...
    [ ...
    0 ,  0 ,  1; ...
    0 ,  0 , -1; ...
    0 ,  1 ,  0 ], ...
    zeros(3), ...
    [ ...
    0 ,  0 , -2; ...
    0 ,  0 ,  0; ...
    0 ,  0 ,  0 ], ...
    zeros(3)];

% tensor hat_G (3xn_el_dof), such that:
% delta_hat_theta = hat_G*shp(hat_Rt)*delta_tilde_a, equation (41)
% equation (43)
hat_G=hat_Csi*hat_Y;

% tensor hat_A (n_el_dofx3), stack(-hat_A_i,I),
% with hat_A_i=spin(V_i+hat_u_i-G)
% stack of spins of position vectors
% equation (52)
hat_A=zeros(n_el_dof,3);
% nodal coordinates of hat triangle (i.e., side-aligned deformed triangle) 
% in element reference frame
coor_hat=coor+hat_u_i;
for i=1:n_el_n
    ind_transl_i=ind_transl(:,i);
    ind_rot_i=ind_rot(:,i);
    hat_A(ind_transl_i,:)=-spin(coor_hat(:,i));
    hat_A(ind_rot_i,:)=eye(3);
end

% tensor L (n_el_dofx3)
% L = stack[I;0]
% equation (52)
L=zeros(n_el_dof,3);
for i=1:n_el_n
    ind_transl_i=ind_transl(:,i);
    L(ind_transl_i,:)=eye(3);
end

% projector hat_P (n_el_dofxn_el_dof), such that:
% delta_hat_a = hat_P*shp(hat_Rt)*delta_tilde_a, equation (50)
% equation (51)
hat_P=eye(n_el_dof)-(L*T+hat_A*hat_G);

end
