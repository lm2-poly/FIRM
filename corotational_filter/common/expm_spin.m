% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function R=expm_spin(theta)

% This function computes the rotation tensor associated with the rotation vector theta
% i.e., it computes R=expm(spin(theta))
% its inverse is the function ax_log
%
% input: 
% 
% theta(3x1)                        rotation vector
%
% output: 
% 
% R(3x3)                            rotation tensor, defined as expm(spin(theta))
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

% adapted from:
% B. Nour-Omid and C.C. Rankin
% Finite rotation analysis and consistent linearization using projectors
% Computer Methods in Applied Mechanics and Engineering 93 (1991) 353-384
% doi: http://dx.doi.org/10.1016/0045-7825(91)90248-5
% and
% C.A. Felippa,  B. Haugen
% A unified formulation of small-strain corotational finite elements: I. Theory
% Comput. Methods Appl. Mech. Engrg. 194 (2005) 2285-2335
% doi: 10.1016/j.cma.2004.07.035

% rotation vector
norm_theta=norm(theta);
spin_theta=spin(theta);
spin_theta_sq=spin_theta*spin_theta;
c1=evaluate('sin(x)/x',norm_theta);
c2=evaluate('(1-cos(x))/x^2',norm_theta);
% Nour-Omid, Rankin, eq. (A.1); Felippa, Haugen, eq. (89)
% equation (1)
R=eye(3)+c1*spin_theta+c2*spin_theta_sq;

end
