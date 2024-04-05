% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [Db, Ds]=linear_isotropic_material(E, nu)

% This routine returns:
% Db(3,3): in-plane reduced (plane stress) stiffness matrix
% Ds(2,2): transversal shear stiffness matrix 
% of an isotropic material with Young modulus E and Poisson's ratio nu
%
% (C) 2010-2013 Paolo Bisegna and Federica Caselli. License: GNU General Public License (GPLv3)


% in-plane reduced (plane stress) stiffness matrix
Db=zeros(3,3);
t1 = (E / (1 + nu)) / 0.2e1;
t2 = E / (1 - nu ^ 2);
t3 = t2 * nu;
Db(1,1) = t2;
Db(1,2) = t3;
Db(2,1) = t3;
Db(2,2) = t2;
Db(3,3) = t1;

% transversal shear stiffness matrix
Ds=zeros(2,2);
Ds(1,1) = t1;
Ds(2,2) = t1;

end
