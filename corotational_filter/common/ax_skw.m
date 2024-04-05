% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function w=ax_skw(W)

% This function computes the axis of the skew-symmetric part of the tensor W
% 
% input
% 
% W(3x3)                    tensor W
% 
% output
% 
% w(3x1)                    axis of the skew-symmetric part of W
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

w=[ W(3,2)-W(2,3); ...
    W(1,3)-W(3,1); ...
    W(2,1)-W(1,2)]/2;

end

