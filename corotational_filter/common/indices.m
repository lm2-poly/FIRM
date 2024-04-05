% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [ind_transl, ind_rot]=indices(n)

% This function returns the indices of translational and rotational dofs
% they are arranged such that:
% a(ind_transl) yields [[u1; v1; w1], ..., [un; vn; wn]],
% i.e., the displacement components (u,v,w) of each node 1,...,n
% and
% a(ind_rot) yields [[q1; r1; s1], ..., [qn; rn; sn]],
% i.e., the rotation vector (q,r,s) of each node 1,...,n
%
% output
%
% ind_transl(3xn)           indices of translational dofs
%
% ind_rot(3xn)              indices of rotational dofs
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

ind=reshape(1:6*n,6,n);
ind_transl=ind(1:3,:);
ind_rot=ind(4:6,:);

end
