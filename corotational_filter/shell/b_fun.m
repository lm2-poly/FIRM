% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [b]=b_fun(coor, bar_u_i, area_def)

% This function computes the vector:
% (1/area_def)*d_(area_def)_d_bar_a
%
% input:
%
% coor(3x3)                 nodal coordinates of reference triangle
%                           coor(:,i) contains the coordinates of node i
%
% bar_u_i(3x3)              core nodal displacements
%                           obtained after filtering out also the rigid rotation chk_R
%                           bar_u_i(:,i) are relevant to node i
%
% area_def                  area of defomred triangle
%
% output:
% 
% b(1x18)                   (1/area_def)*d_(area_def)_d_bar_a
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

% computation based on nodal coordinates of the bar-triangle
% in the frame (G,e,h,n)
bar_coor=coor+bar_u_i;
bar_x=bar_coor(1,:);
bar_y=bar_coor(2,:);
% equation (106)
b=[ ...
    bar_y(2)-bar_y(3), bar_x(3)-bar_x(2), 0, ...
    0, 0, 0, ...
    bar_y(3)-bar_y(1), bar_x(1)-bar_x(3), 0, ...
    0, 0, 0, ...
    bar_y(1)-bar_y(2), bar_x(2)-bar_x(1), 0, ...
    0, 0, 0, ...
    ]/(2*area_def);

end
