% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [alpha_b, beta, beta0]=OPT_parameters(nu)

% This functions computes the 11 optimal parameters of OPT membrane for an isotropic material
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

% Optimal ANDES Template

% optimal alpha_b and beta, eq. (42)
alpha_b=3/2;
beta=[1,2,1,0,1,-1,-1,-1,-2];

% equation (42)
beta0=(1-4*nu^2)/2;

% beta0 must not go below a positive treshold, say 0.01
% see text following eq. (42)
beta0=max(beta0,0.01);

end
