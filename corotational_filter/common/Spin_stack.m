% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [V]=Spin_stack(v,n)

% This function restuns the stack of skew-symmetrix tensors 
% associated with the 2n 3-component vectors contained in v
%
% input:
% v(6nx1)	stack of 2n 3-component vectors
%
% output:
% V(6nx3)   stack of the skew-symmetrix tensors associated with the six 3-component vectors contained in v
%           i.e., V=[spin(v(1:3)); spin(v(4:6)); ... ; spin(v(6n-2:6n))]
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

% initialize V
V=zeros(6*n,3);
% initialize pointer
ind=0;
% equation (86)
% cycle on the 2n 3-component vectors
for i=1:2*n
    % convert the 3-component vector v(range) to skew-symmetrix tensor 
    range=ind+1:ind+3;
    V(range,:)=spin(v(range));
    % update pointer
    ind=ind+3;
end

end
