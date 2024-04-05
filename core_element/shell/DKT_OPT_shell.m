% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [ ... % material tangent stiffness matrix
    K_t_mat, ...
    ... % nodal forces equivalent to internal forces
    q_int, ...
    ... % opposite of nodal forces equivalent to applied loads
    q_load, Q_load, QM_load, H_load, ...
    ... % resultants of applied loads
    f_load, m_G_load, M_G_load]= DKT_OPT_shell( ...
    ... % nodal coordinates and nodal parameters
    coor, a, ...
    ... % number of load conditions and applied loads
    n_load, load, ...
    ... % material parameters
    mat)

% This function implements the triangular shell element: OPT/CST membrane + DKT plate
% linear elastic material
%
% tiangle located in the x,y plane
%
% shell element dofs:
%   ux, uy, uz, rx, ry, rz at node 1, then at node 2, then at node 3
%
% input:
%
% coor(2x3)                 nodal coordinates
%       coor(:,i):          (x,y) coordinates of node i
%
% a(18x1)                   nodal parameters
%                           a=(u_1;r_1;u_2;r_2;u_3;r_3)
%                           u_i=(ux,uy,uz) displacements of node i
%                           r_i=(rx,ry,rz) rotations of node i
%
% n_load                    number of load conditions
%
% load(6x3xn_load)          distributed loads per unit area, linear over the triangle
%       load(1:3,i,lc):     distributed forces (lx,ly,lz) at node i, load condition lc
%       load(4:6,i,lc):     distributed couples (mx,my,mz) at node i, load condition lc
%
% mat:                      geometric and material parameters; structure with the following fields:
%
%     n_g                   number of Gauss points used for numeric integration over triangle
%                           admissible values are: 1, 3, 4, 7, 12
%     constitutive          constitutive law: struct with the following fields
%         type              type of constitutive law, string
%                           presently available:
%                           'linear_elastic_isotropic_plane_stress' alias 'saint_venant_kirchhoff_isotropic_plane_stress'
%         param             parameters of constitutive law,
%                           for 'linear_elastic_isotropic_plane_stress' alias 'saint_venant_kirchhoff_isotropic_plane_stress':
%                           E, nu
%     th                    thickness
%
% output:
%
% K_t_mat(18x18)            (small displacement) (linear) material tangent stiffness matrix
%
% q_int(18x1)               nodal forces equivalent to internal forces
%
% q_load(18xn_load)         opposite of nodal forces equivalent to applied loads
%                           -int( (displacement shape functions)'*(distributed forces) + (rotation shape functions)'*(distributed couples))
%
% Q_load(18x3xn_load)       -int( (displacement shape functions)'*spin(distributed forces))
%
% QM_load(18x3xn_load)      -int( (displacement shape functions)'*spin(distributed forces) + (rotation shape functions)'*spin(distributed couples) )
%
% H_load(18x18xn_load)      -int( Delta_(displacement shape functions)'*(distributed forces) + Delta_(rotation shape functions)'*(distributed couples) )
%
% f_load(3xn_load)          resultant force: int( distributed forces )
%
% m_G_load(3xn_load)        resultant moment wrt G: int( spin(deformed position vector wrt G)*(distributed forces) + (distributed couples) )
%
% M_G_load(3x3xn_load)      int( spin(deformed position vector wrt G)*spin(distributed forces) + spin(distributed couples) )
%
% numeric integration is implemented here for clarity;
% however, closed-form expressions of stiffness matrices and load contributions are easily obtained
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

% extract geometric and material parameters
n_g          = mat.n_g; % number of Gauss points used for numeric integration over triangle
th           = mat.th; % tickness
constitutive = mat.constitutive; % type of constitutive law and relevant parameters

% extract type of constitutive law and params
type        = constitutive.type;
param       = constitutive.param;

switch type
    case {'linear_elastic_isotropic_plane_stress', 'saint_venant_kirchhoff_isotropic_plane_stress'}
        % extract material parameters
        E            = param.E;
        nu           = param.nu;
        
        % Db(3,3): in-plane reduced (plane stress) stiffness matrix
        [Db, ~]=linear_isotropic_material(E, nu);
        
    otherwise
        error(['constitutive law ',type, ' not allowed in linear shell. Try material nonlinear shell'])
end

% membranal stiffness
AA = th*Db;
% bending stiffness
DD = th^3/12*Db;

% nodal coordinates
x=coor(1,:);
y=coor(2,:);

% centroid
G=sum([x; y],2)/3;

% initialize element matrices
% (small displacement) (linear) material tangent stiffness matrix
K_t_mat=zeros(18);

% nodal forces equivalent to distributed loads
% -int( (displacement shape functions)'*(distributed forces) + (rotation shape functions)'*(distributed couples))
q_load=zeros(18,n_load);
% -int( (displacement shape functions)'*spin(distributed forces) )
Q_load=zeros(18,3,n_load);
% -int( (displacement shape functions)'*spin(distributed forces) + (rotation shape functions)'*spin(distributed couples) )
QM_load=zeros(18,3,n_load);
% -int( Delta_(displacement shape functions)'*(distributed forces) + Delta_(rotation shape functions)'*(distributed couples) )
H_load=zeros(18,18,n_load);  % this is zero, since shape functions are linear in nodal dofs

% resultant force: int( distributed forces )
f_load=zeros(3,n_load);
% resultant moment wrt G: int( spin(deformed position vector wrt G)*(distributed forces) + (distributed couples) )
m_G_load=zeros(3,n_load);
% int( spin(deformed position vector wrt G)*spin(distributed forces) + spin(distributed couples) )
M_G_load=zeros(3,3,n_load);


% 9 membrane dofs
% ux, uy, psi at node 1, then
% ux, uy, psi at node 2, then
% ux, uy, psi at node 3
% mapping to 18 shell dofs
membr_dof=repmat([true; true; false; false; false; true],3,1);

% shape function of membrane

% M_etax(3,9), M_etay(3,9), M_gamxy(3,9), M_ux(6,9), M_uy(6,9), M_psi(3,9)

% the shape functions of in-plane deformations etax, etay and gamxy
% are linear homogeneous polynomials of areal coordinates L=[1-xi-eta, xi, eta];
% N_etax (shape function of etax) is given by L*M_etax,
% analogously, N_etay=L*M_etay, N_gamxy=L*M_gamxy;

% the shape function of in-plane displacements ux and uy are quadratic homogeneous polynomial of areal coordinates
% N_ux (shape function of ux) is given by Lsq*M_ux;
% analogously, N_uy=Lsq*M_uy

% the shape function of drilling rotation psi
% is a linear homogeneous polynomial of areal coordinates
% N_psi (shape function of psi) is given by L*M_psi

% OPT membrane:
% adapted from:
% Carlos A. Felippa, "A study of optimal membrane triangles with drilling freedoms"
% Comput. Methods Appl. Mech. Engrg. 192(2003) 2125-2168
% dx.doi.org/10.1016/S0045-7825

% shape functions of OPT membrane
[M_etax, M_etay, M_gamxy, M_ux, M_uy, M_psi]=OPT_shape_functions(x, y, nu);

% 9 plate dofs
% w, theta_x, theta_y at node 1, then
% w, theta_x, theta_y at node 2, then
% w, theta_x, theta_y at node 3
% mapping to 18 shell dofs
plate_dof=repmat([false; false; true; true; true; false],3,1);

% DKT plate

% shape functions of DKT plate
%
% M_thetax(6,9); M_thetay(6,9); M_w(15,9)

% the shape functions of thetax and thetay
% are quadratic homogeneous polynomials of areal coordinates L=[1-xi-eta, xi, eta];
% N_thetax (shape function of thetax) is obtained by Lsq*M_thetax
% where Lsq contains the quadratic homogeneous monomials of areal coordinates
% analogously, N_thetay=Lsq*M_thetay

% the shape function of w is a quartic homogeneous polynomial of areal coordinates
% N_w (shape function of w) is obtained by Lsqsq*M_w, where:
% where Lsqsq contains the quartic homogeneous monomials of areal coordinates

[M_thetax, M_thetay, M_w]=DKT_shape_functions(x, y);


% membrane 9x9 stiffness matrix
K_membr=zeros(9);
% plate 9x9 stiffness matrix
K_plate=zeros(9);

% all-zero vector
zero=zeros(1,9);

% Gaussian quadrature for a triangular domain
[L1, L2, L3, w]=gauss_tri(n_g);

% linear map (jacobian determinant and jacobian inverse)
[~, detJ, Jinvt]=lin_tri_map(x, y);

% cycle on Gauss points
for g=1:n_g
    
    % area element
    w_detJ=w(g)*detJ/2;
    
    %  areal coordinates
    L=[L1(g), L2(g), L3(g)];
    
    % quadratic monomials of areal coordinates
    Lsq=[L(1)^2, L(2)^2, L(3)^2, L(2)*L(3), L(3)*L(1), L(1)*L(2)];
    
    % derivatives (wrt xi,eta) of quadratic monomials of areal coordinates
    d_Lsq_dxi=[ ...
        [-2*L(1), 2*L(2), 0,      L(3), -L(3),     L(1)-L(2)]; ...
        [-2*L(1), 0,      2*L(3), L(2), L(1)-L(3), -L(2)    ]; ...
        ];
    
    % quartic monomials of areal coordinates
    Lsqsq = [ ...
        L(1) ^ 4, ...
        L(2) ^ 4, ...
        L(3) ^ 4, ...
        L(2) ^ 3 * L(3), ...
        L(2) * L(3) ^ 3, ...
        L(3) ^ 3 * L(1), ...
        L(3) * L(1) ^ 3, ...
        L(1) ^ 3 * L(2), ...
        L(1) * L(2) ^ 3, ...
        L(2) ^ 2 * L(3) ^ 2, ...
        L(3) ^ 2 * L(1) ^ 2, ...
        L(1) ^ 2 * L(2) ^ 2, ...
        L(1) ^ 2 * L(2) * L(3), ...
        L(2) ^ 2 * L(3) * L(1), ...
        L(3) ^ 2 * L(1) * L(2)];
    
    
    % membrane shape functions
    
    % shape functions of in-plane displacements and drilling rotation
    N_u=[Lsq*M_ux; Lsq*M_uy];
    N_psi=L*M_psi;
    
    % shape functions of in-plane deformations etax, etay, gamxy
    N_eta=[L*M_etax; L*M_etay; L*M_gamxy];
    
    % plate shape functions
    
    % shape functions of bending rotations thetax and thetay
    N_theta=[Lsq*M_thetax; Lsq*M_thetay];
    
    % shape functions of derivative of thetax wrt xi,eta
    d_N_thetax_dxi=d_Lsq_dxi*M_thetax;
    % shape functions of derivative of thetax wrt x,y,
    d_N_thetax_dx=Jinvt*d_N_thetax_dxi;
    
    % shape functions of derivative of thetay wrt xi,eta
    d_N_thetay_dxi=d_Lsq_dxi*M_thetay;
    % shape functions of derivative of thetay wrt x,y
    d_N_thetay_dx=Jinvt*d_N_thetay_dxi;
    
    % shape functions of curvatures kappax, kappay, kappaxy
    N_kappa=zeros(3,9);
    N_kappa(1,:)= d_N_thetay_dx(1,:);
    N_kappa(2,:)=-d_N_thetax_dx(2,:);
    N_kappa(3,:)= d_N_thetay_dx(2,:)-d_N_thetax_dx(1,:);
    
    % shape functions of w
    N_w=Lsqsq*M_w;
    
    % stiffness matrix
    
    % membrane contribution
    % int_thickness 1/2*(eta*Db*eta)
    K_membr=K_membr+(N_eta'*AA*N_eta)*w_detJ;
    % plate contribution
    % int_thickness 1/2*z^2*(kappa*Db*kappa)
    K_plate=K_plate+(N_kappa'*DD*N_kappa)*w_detJ;
    
    
    % distributed loads per unit area
    
    % position vector wrt centroid G in the reference configuration
    rxy=[x; y]*L'-G;
    rz=0;
    
    % position vector wrt centroid G in the deformed configuration
    rxy_def=rxy+N_u*a(membr_dof);
    rz_def=rz+N_w*a(plate_dof);
    r_def=[rxy_def; rz_def];
    spin_r_def=spin(r_def);
    
    % transposed matrices of shape functions
    N_u_zero_t=[N_u; zero]';
    N_theta_zero_t=[N_theta; zero]';
    zero_zero_N_psi_t=[zero; zero; N_psi]';
    zero_zero_N_w_t=[zero; zero; N_w]';
    
    % cycle on load types
    for j=1:n_load
        % load (force and couple) intensity at Gauss point, times Gauss weight
        l_w_detJ=load(1:3,:,j)*L'*w_detJ;
        m_w_detJ=load(4:6,:,j)*L'*w_detJ;
        
        % computes load contribution (if loads are present)
        if any(l_w_detJ) || any(m_w_detJ)
            
            % spin of load vectors
            spin_l_w_detJ=spin(l_w_detJ);
            spin_m_w_detJ=spin(m_w_detJ);
            
            % opposite of nodal forces equivalent to applied loads, i.e.,
            % -int( (displacement shape functions)'*(distributed forces) + (rotation shape functions)'*(distributed couples) )
            q_load(membr_dof,j)=q_load(membr_dof,j) -(N_u_zero_t*l_w_detJ+zero_zero_N_psi_t*m_w_detJ);
            q_load(plate_dof,j)=q_load(plate_dof,j) -(zero_zero_N_w_t*l_w_detJ+N_theta_zero_t*m_w_detJ);
            
            % -int( (displacement shape functions)'*spin(distributed forces) )
            Q_load(membr_dof,:,j)=Q_load(membr_dof,:,j) -N_u_zero_t*spin_l_w_detJ;
            Q_load(plate_dof,:,j)=Q_load(plate_dof,:,j) -zero_zero_N_w_t*spin_l_w_detJ;
            
            % -int( (displacement shape functions)'*spin(distributed forces) + (rotation shape functions)'*spin(distributed couples) )
            QM_load(membr_dof,:,j)=QM_load(membr_dof,:,j) -(N_u_zero_t*spin_l_w_detJ+zero_zero_N_psi_t*spin_m_w_detJ);
            QM_load(plate_dof,:,j)=QM_load(plate_dof,:,j) -(zero_zero_N_w_t*spin_l_w_detJ+N_theta_zero_t*spin_m_w_detJ);
            
            % load resultants
            
            % resultant force: int( (distributed forces) )
            f_load(:,j)=f_load(:,j) + l_w_detJ;
            
            % resultant moment wrt G: int( spin(position vector wrt G)*(distributed forces) + (distributed couples) )
            m_G_load(:,j)=m_G_load(:,j) + (spin_r_def*l_w_detJ+m_w_detJ);
            
            % int( spin(position vector wrt G)*spin(distributed forces) +  spin(distributed couples) )
            M_G_load(:,:,j)=M_G_load(:,:,j) +(spin_r_def*spin_l_w_detJ+spin_m_w_detJ);
            
        end
        
    end
    
end

% assemble into 18x18 shell matrices

% stiffness matrix
K_t_mat(membr_dof,membr_dof)=K_membr;
K_t_mat(plate_dof,plate_dof)=K_plate;

% nodal forces equivalent to internal forces
q_int=K_t_mat*a;

end
