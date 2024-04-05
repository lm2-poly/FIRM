% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [shp_R]=shp(R, n)

% This function returns the block diagonal tensor collecting 6 copies of R
% i.e., shp_R=blkdiag(R,...,R) (2n copies of R)
%
% input
% 
% R(3x3)                    rotation tensor
%
% n                         number of nodes
%
% output
%
% shp_R(6nx6n)              block-diagonal assembly of 2n copies of R
%
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

% initialize shp_R
shp_R=zeros(6*n);
% initialize pointer
ind=0;
for i=1:2*n
    range=ind+1:ind+3;
    % copy R on the diagonal of shp_R
    shp_R(range,range)=R;
    % update pointer
    ind=ind+3;
end

end
