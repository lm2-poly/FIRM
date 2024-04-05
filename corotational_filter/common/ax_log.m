% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function theta=ax_log(R)

% This function computes the axis of the log of the rotation tensor R
% i.e., it computes theta such that R=expm(spin(theta))
% its inverse is the function expm_spin
%
% input:
%
% R(3x3)                    rotation tensor
%
% output:
%
% theta(3x1)                rotation vector, such that R=expm(spin(theta));
%
% robust algorithm in:
% Comment on "Singularity-Free Extraction of a Quaternion from a Direction-Cosine Matrix"
% Richard A. Spurrier,
% JOURNAL OF SPACECRAFT AND ROCKETS 1978
% 0022-4650 vol.15 no.4 (255-255)
% doi: 10.2514/3.57311
% cf. also
% C.A. Felippa,  B. Haugen
% A unified formulation of small-strain corotational finite elements: I. Theory
% Comput. Methods Appl. Mech. Engrg. 194 (2005) 2285-2335
% doi: 10.1016/j.cma.2004.07.035, eq. (97)
%
% the naive procedure:
% [psi]=ax_skw(R);
% norm_psi=norm(psi);
% c=evaluate('asin(x)/x',norm_psi);
% theta=c*psi;
% fails when abs(norm(theta)) is greater that pi/2
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

% initialize quaternions
p=zeros(3,1);

% trace of R
% trace(R) is slower
tr_R=R(1,1)+R(2,2)+R(3,3);

% algebraically largest diagonal term
[R_ii,i_max]=max(diag(R));
% compare tr_R with diagonal terms of R
if tr_R>R_ii
    % tr_R is algebraically larger than R_ii
    p0=sqrt(1+tr_R)/2;
    for i=1:3
        % construct i,j,k as cyclic permutation of 1,2,3
        j=mod(i,3)+1;
        k=mod(j,3)+1;
        % compute quaternions
        p(i)=(R(k,j)-R(j,k))/(4*p0);
    end
else
    % R_ii is algebraically larger than tr_R
    i=i_max;
    j=mod(i,3)+1;
    k=mod(j,3)+1;
    % compute quaternions
    p(i)=sqrt(R_ii/2+(1-tr_R)/4);
    four_p_i=4*p(i);
    p(j)=(R(j,i)+R(i,j))/four_p_i;
    p(k)=(R(k,i)+R(i,k))/four_p_i;
    p0=(R(k,j)-R(j,k))/four_p_i;
    % in case p0 turns out to be negative, change sign, in order to assure
    % that norm(theta) will not be greater than pi
    if p0<0
        p0=-p0;
        p=-p;
    end
end

% compute theta from quaternions
% norm of p
norm_p=norm(p);
% sin(norm(theta)/2), given by norm_p
sin_norm_th2=norm_p;
% cos(norm(theta)/2), given by p0
cos_norm_th2=p0;
% norm(theta)/2, via atan2
norm_th2=atan2(sin_norm_th2,cos_norm_th2);

% theta, proportional to p, with length norm_theta
% when R=eye(3), the naive formula: theta=2*norm_th2/norm_p*p returns nans;
% use Taylor expansion, instead
theta=2/evaluate('sin(x)/x',norm_th2)*p;

end
