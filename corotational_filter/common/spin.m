% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function W=spin(w)

% This function computes the skew-symmetrix tensor associated with a 3-component vector
%
% input:
% w(3x1)    3-component vector
%
% output:
% W(3x3)    skew-symmetrix tensor associated with the vector w
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

% elements of W, arranged columnwise
W=[0, w(3), -w(2), -w(3), 0, w(1), w(2), -w(1), 0];
% reshape into a 3x3 skew tensor
W=reshape(W,3,3);

end
