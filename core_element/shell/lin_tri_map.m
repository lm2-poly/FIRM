% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [J,detJ,Jinvt]=lin_tri_map(x,y)

% This routine computes the linear map from parent triangle to physical triangle:
% it returns jacobian, jacobian determinant and jacobian inverse-transpose
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

% jacobian
J=[ x(2)-x(1), x(3)-x(1) ; ...
    y(2)-y(1), y(3)-y(1)];

% jacobian determinant
% detJ = -x(2) * y(1) + x(2) * y(3) - x(1) * y(3) + x(1) * y(2) - x(3) * y(2) + x(3) * y(1);
detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);

% inverse transpose jacobian
Jinvt=zeros(2,2);
t9 = 0.1e1 / (detJ);

% Jinvt(1,1) = -(y(1) - y(3)) * t9;
% Jinvt(1,2) = (-y(2) + y(1)) * t9;
% Jinvt(2,1) = (x(1) - x(3)) * t9;
% Jinvt(2,2) = -(-x(2) + x(1)) * t9;

Jinvt(1,1) =  J(2,2) * t9;
Jinvt(1,2) = -J(2,1) * t9;
Jinvt(2,1) = -J(1,2) * t9;
Jinvt(2,2) =  J(1,1) * t9;

end
