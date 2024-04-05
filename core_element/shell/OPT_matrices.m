% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [B, T_theta_u, Te, Q1, Q2, Q3, A]=OPT_matrices(alpha_b, beta, x, y)

% This functions computes useful matrices for OPT membrane
%
% adapted from:
% Carlos A. Felippa, "A study of optimal membrane triangles with drilling freedoms"
% Comput. Methods Appl. Mech. Engrg. 192(2003) 2125-2168
% dx.doi.org/10.1016/S0045-7825
% Equation numbers below refer to that paper
%
% B(3x9)=L'/V, where L is the force-lumping matrix, eq. (22)
% T_theta_u(3x9): hierarchical rotations, eq. (16)
% Te(3x3): 'straingage rosette' trasformation, eq. (12)
% Qi(3x3), i=1,2,3: matrix relating the natural strains eps_i at corner i
%   to the deviatoric corner curvatures theta_tilda, eq. (29)
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

% side projections along the coordinate axes
% x
x12=x(1)-x(2);
x23=x(2)-x(3);
x31=x(3)-x(1);
x21=-x12;
x32=-x23;
x13=-x31;
% y
y12=y(1)-y(2);
y23=y(2)-y(3);
y31=y(3)-y(1);
y21=-y12;
y32=-y23;
y13=-y31;

% element area 2*A=norm(cross(l12,l13))
A=(y21*x13-x21*y13)/2;

% B(3x9), equal to L'/V
% L(9x3) is the force-lumping matrix, eq. (22)
% given a constant stress vector s={sxx,syy,sxy}, L*s is a 9-component
% vector containing the resulting forces and couples at the element nodes
% on the other hand, B is a 3x9 matrix which can be interpreted as a
% strain-displacement matrix: indeed, B times nodal displacement {d}
% gives the constant mean strain inside the elements indiced by {d};
% in the case alpha_b=0, B reduces to the CST-3/6C matrix
B= [[y23;           0;              x32], ...
    [0;             x32;            y23], ...
    [y23*(y13-y21); x32*(x31-x12); (x31*y13-x12*y21)*2]*alpha_b/6, ...
    [y31;           0;              x13], ...
    [0;             x13;            y31], ...
    [y31*(y21-y32); x13*(x12-x23); (x12*y21-x23*y32)*2]*alpha_b/6, ...
    [y12;           0;              x21], ...
    [0;             x21;            y12], ...
    [y12*(y32-y13); x21*(x23-x31); (x23*y32-x31*y13)*2]*alpha_b/6]/(2*A);

% T_theta_u(3x9): hierarchical rotations, eq. (16)
% T_theta_u is used to extract the hierarchical corner rotations
% from the total corner rotation, subtract the mean or CST
% (constant strain triangle) rotation
T_theta_u=1/(4*A)*[x32,y32,0,x13,y13,0,x21,y21,0];
T_theta_u=[T_theta_u;T_theta_u;T_theta_u];
T_theta_u(:,[3,6,9])=eye(3);

% squared sides lenghts: LLij=xij^2 + yij^2), LLij=LLji
LL21=x21^2+y21^2;
LL32=x32^2+y32^2;
LL13=x13^2+y13^2;

% Te(3x3): 'straingage rosette' trasformation, eq. (12)
% These are extensional (direct) strains along the three side
% directions intrinsically related to the triangle geometry
% The natural strain are related to Cartesian strains by the
% 'straingage rosette' trasformation:
% [eps12; eps32; eps13]= inv(Te)* [exx; eyy; 2exy]
% Te is constant over the triangle
Te=[[y23*y13; x23*x13; (y23*x31+x32*y13)]*LL21, ...
    [y31*y21; x31*x21; (y31*x12+x13*y21)]*LL32, ...
    [y12*y32; x12*x32; (y12*x23+x21*y32)]*LL13]/(4*A^2);

% % %         % check, eq. (11)
% % %         Te_inv=repmat([1/LL21;1/LL32;1/LL13],1,3).*[x21^2, y21^2, x21*y21;
% % %             x32^2, y32^2, x32*y32;
% % %             x13^2, y13^2, x13*y13];
% % %         norm(Te-inv(Te_inv))/norm(Te)

% % %         % displacements and hierarchical rotations, eq. (17)
% % %         T_R=eye(9);
% % %         T_R([3,6,9],:)=T_theta_u;
% % %         % check, eq. (23)
% % %         norm(B*T_R-B)/norm(B)

% Qi(3x3), i=1,2,3: matrix relating the natural strains eps_i at corner i
% to the deviatoric corner curvatures theta_tilda, eq. (29)
M_LL=[1/LL21; 1/LL32; 1/LL13]*(2*A/3);
M_LL=[M_LL,M_LL,M_LL];

% first equation in (29)
M_beta1=reshape(beta([1,4,7, 2,5,8, 3,6,9]),3,3);
Q1=M_LL.*M_beta1;

% second equation in (29)
M_beta2=reshape(beta([9,3,6, 7,1,4, 8,2,5]),3,3);
Q2=M_LL.*M_beta2;

% third equation in (29): note that there is a typo in that equation, third row
M_beta3=reshape(beta([5,8,2, 6,9,3, 4,7,1]),3,3);
Q3=M_LL.*M_beta3;

% % % % check: Q1+Q2+Q3=0
% % % Q1+Q2+Q3

end
