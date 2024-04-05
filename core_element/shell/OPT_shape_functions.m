% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [M_etax, M_etay, M_gamxy, M_ux, M_uy, M_psi]=OPT_shape_functions(x, y, nu)

% This functions computes the shape functions of OPT membrane
%
% input: 
% 
% x(3), y(3)                nodal coordinates (the triangle is located in the x,y plane)
%
% nu                        Poisson's ratio
%
% th                        membrane thickness
%
% output:
%
% M_etax(3,9), M_etay(3,9), M_gamxy(3,9)
% are used to compute the shape functions of in-plane deformations etax, etay and gamxy,
% which are linear homogeneous polynomials of areal coordinates L=[1-xi-eta, xi, eta];
% N_etax (shape function of etax) is given by L*M_etax,
% analogously, N_etay=L*M_etay, N_gamxy=L*M_gamxy; 
%
% M_ux(6,9), M_uy(6,9)
% are used to compute the shape function of in-plane displacements ux and uy,
% which are quadratic homogeneous polynomial of areal coordinates;
% N_ux (shape function of ux) is given by Lsq*M_ux, where:
% Lsq=[L(1)^2, L(2)^2, L(3)^2, L(2)*L(3), L(3)*L(1), L(1)*L(2)]
% analogously, N_uy=Lsq*M_uy
%
% M_psi(3,9) 
% is used to compute the shape function of drilling rotation psi,
% which is a linear homogeneous polynomial of areal coordinates
% N_psi (shape function of psi) is given by L*M_psi
%
% the 9 dofs are ordered as follows: 
%   ux, uy, psi at node 1, then at node 2, then at node 3
%
% adapted from:
% Carlos A. Felippa, "A study of optimal membrane triangles with drilling freedoms"
% Comput. Methods Appl. Mech. Engrg. 192(2003) 2125-2168
% dx.doi.org/10.1016/S0045-7825
% Equation numbers below refer to that paper
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

% compute the 11 optimal parameters for an orthotropic material
[alpha_b, beta, beta0]=OPT_parameters(nu);

% compute useful matrices:
% B(3x9)=L'/V, where L is the force-lumping matrix, eq. (22)
% T_theta_u(3,9): hierarchical rotations, eq. (16)
% Te(3,3): 'straingage rosette' trasformation, eq. (12)
% Qi(3,3), i=1,2,3: matrix relating the natural strains eps_i at corner i 
%   to the deviatoric corner curvatures theta_tilda, eq. (29)
% A: element area
% V: element volume
[B, T_theta_u, Te, Q1, Q2, Q3, ~]=OPT_matrices(alpha_b, beta, x, y);

% in-plane deformations etax, etay, gamxy
% the coefficient beta0 in eq. (49) does NOT coincide with the one used in eq. (32):
% indeed, in order to obtain eq. (32) integrating the elastic energy
% computed from the deformations in eq. (49), it is necessary that in the
% last equation beta0 is replaced by beta0_tilda, given by:
beta0_tilda=3/2*sqrt(beta0);

% [etax;etay;gamxy] = M1*L(1)+M2*L(2)+M(3)*L(2), eq. (49)
% the constant part B is reduced to a homogeneous polynomial 
% in the areal coordinates using: 1=L(1)+L(2)+L(3)
% M1(3,9), M2(3,9), M3(3,9) are given by:
M1=B+(beta0_tilda*Te*Q1*T_theta_u);
M2=B+(beta0_tilda*Te*Q2*T_theta_u);
M3=B+(beta0_tilda*Te*Q3*T_theta_u);

% computing matrices M_etax(3,9), M_etay(3,9), M_gamxy(3,9);
% the shape functions of in-plane deformations etax, etay and gamxy
% are linear homogeneous polynomials of areal coordinates L=[1-xi-eta, xi, eta];
% N_etax (shape function of etax) is given by L*M_etax,
% analogously, N_etay=L*M_etay, N_gamxy=L*M_gamxy; 
M_etax =[M1(1,:);M2(1,:);M3(1,:)];
M_etay =[M1(2,:);M2(2,:);M3(2,:)];
M_gamxy=[M1(3,:);M2(3,:);M3(3,:)];


% the displacement field of OPT membrane is not defined
% an Allman-type displacement field is used instead;
% the shape function of in-plane displacements ux and uy are quadratic homogeneous polynomial of areal coordinates
% N_ux (shape function of ux) is given by Lsq*M_ux, where:
% Lsq=[L(1)^2, L(2)^2, L(3)^2, L(2)*L(3), L(3)*L(1), L(1)*L(2)]
% analogously, N_uy=Lsq*M_uy
M_ux = [ ...
    1 0 0 0 0 0 0 0 0; ...
    0 0 0 1 0 0 0 0 0; ...
    0 0 0 0 0 0 1 0 0; ...
    0 0 0 1 0 alpha_b * (-y(3) + y(2)) / 2 1 0 -alpha_b * (-y(3) + y(2)) / 2; ...
    1 0 alpha_b * (y(1) - y(3)) / 2 0 0 0 1 0 -alpha_b * (y(1) - y(3)) / 2; ...
    1 0 alpha_b * (-y(2) + y(1)) / 2 1 0 -alpha_b * (-y(2) + y(1)) / 2 0 0 0];

M_uy = [ ...
    0 1 0 0 0 0 0 0 0; ...
    0 0 0 0 1 0 0 0 0; ...
    0 0 0 0 0 0 0 1 0; ...
    0 0 0 0 1 -alpha_b * (-x(3) + x(2)) / 2 0 1 alpha_b * (-x(3) + x(2)) / 2; ...
    0 1 -alpha_b * (x(1) - x(3)) / 2 0 0 0 0 1 alpha_b * (x(1) - x(3)) / 2; ...
    0 1 -alpha_b * (-x(2) + x(1)) / 2 0 1 alpha_b * (-x(2) + x(1)) / 2 0 0 0];

% the shape function of drilling rotation psi
% is assumed to be a linear homogeneous polynomial of areal coordinates
% N_psi (shape function of psi) is given by L*M_psi
M_psi=[ ...
    0 0 1 0 0 0 0 0 0; ...
    0 0 0 0 0 1 0 0 0; ...
    0 0 0 0 0 0 0 0 1];

end
