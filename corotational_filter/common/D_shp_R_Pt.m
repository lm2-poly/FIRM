% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [K_shp_R_Pt]=D_shp_R_Pt(v, shp_R, G, P, flt_I, n)

% This function computes the tensor:
% d_(shp_R*P'*v)_d_a, with constant v
%
% input:
%
% v(6nx1)                   constant vector
%
% shp_R(6nx6n)              block-diagonal assembly of 2n copies of R
%
% G(3x6n), P(6nx6n)         tensors needed to filter out the rigid rotation R
%
% flt_I(6nx6n)              block-diagonal assembly (I,0,...) n times
%
% output:
%
% K_shp_R_Pt(6nx6n)         tensor d_(shp_R*P'*v)_d_a, with constant v
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

% compute Spin_stack(P'*v, n)*G
Spin_stack_Pt_v_G=Spin_stack(P'*v, n)*G;

% compute d_(shp_R*P'*v)_d_a, with constant v
% equation (125) or (126)
K_shp_R_Pt = -shp_R*( ...
    ... % contribution from derivative of shp_R
    + Spin_stack_Pt_v_G ...
    ... % contribution from derivative of P'
    + Spin_stack_Pt_v_G'*flt_I*P ...
    )*shp_R';

end
