% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [L1, L2, L3, w]=gauss_tri(n_g)

% This routine returns the Gauss points and weights for parent triangle
%
% presently implemented n_g values:
%	 1:       1-point gauss integration
%	 3:       3-point integration: mid-edge points
%	 4:       4-point gauss integration
%	 7:       7-point gauss integration
%	12:      12-point order 6 formula

% initialize
g=zeros(4,n_g);

% useful constants
ww=0.3333333333333333;

switch n_g
    
    case 1
        % 1-point gauss integration
        
        g(1,1) = ww;
        g(2,1) = ww;
        g(3,1) = ww;
        g(4,1) = 1.;
        
    case 3
        % 3-point integration: mid-edge points
        
        g(1,1) = 0.;
        g(2,1) = 0.5;
        g(3,1) = 0.5;
        g(4,1) = ww;
        
        g(1,2) = 0.5;
        g(2,2) = 0.;
        g(3,2) = 0.5;
        g(4,2) = ww;
        
        g(1,3) = 0.5;
        g(2,3) = 0.5;
        g(3,3) = 0.;
        g(4,3) = ww;
        
    case 4
        % 4-point gauss integration
        
        g(1,1) =  ww;
        g(2,1) =  ww;
        g(3,1) =  ww;
        g(4,1) = -27./48.;
        
        g(1,2) =  0.6;
        g(2,2) =  0.2;
        g(3,2) =  0.2;
        g(4,2) =  25./48.;
        
        g(1,3) =  0.2;
        g(2,3) =  0.6;
        g(3,3) =  0.2;
        g(4,3) =  g(4,2);
        
        g(1,4) =  0.2;
        g(2,4) =  0.2;
        g(3,4) =  0.6;
        g(4,4) =  g(4,2);
        
    case 7
        % 7-point gauss integration
        
        r0      =  sqrt(15.0);
        r1      =  3./7.;
        r2      =  (r0 + r0)/21.;
        
        g(1,1) =  ww;
        g(2,1) =  g(1,1);
        g(3,1) =  g(1,1);
        g(4,1) =  0.225;
        
        g(1,2) =  r1 + r2;
        g(2,2) =  0.5 - 0.5*g(1,2);
        g(3,2) =  g(2,2);
        g(4,2) =  (155. - r0)/1200.;
        
        g(1,3) =  g(2,2);
        g(2,3) =  g(1,2);
        g(3,3) =  g(2,2);
        g(4,3) =  g(4,2);
        
        g(1,4) =  g(2,2);
        g(2,4) =  g(2,2);
        g(3,4) =  g(1,2);
        g(4,4) =  g(4,2);
        
        g(1,5) =  r1 - r2;
        g(2,5) =  0.5 - 0.5*g(1,5);
        g(3,5) =  g(2,5);
        g(4,5) =  (155. + r0)/1200.;
        
        g(1,6) =  g(2,5);
        g(2,6) =  g(1,5);
        g(3,6) =  g(2,5);
        g(4,6) =  g(4,5);
        
        g(1,7) =  g(2,5);
        g(2,7) =  g(2,5);
        g(3,7) =  g(1,5);
        g(4,7) =  g(4,5);
        
    case 12
        % 12-point order 6 formula
        
        g(1, 1) = 0.873821971016996;
        g(2, 1) = 0.063089014491502;
        g(3, 1) = 0.063089014491502;
        g(4, 1) = 0.050844906370207;
        
        g(1, 2) = 0.063089014491502;
        g(2, 2) = 0.873821971016996;
        g(3, 2) = 0.063089014491502;
        g(4, 2) = 0.050844906370207;
        
        g(1, 3) = 0.063089014491502;
        g(2, 3) = 0.063089014491502;
        g(3, 3) = 0.873821971016996;
        g(4, 3) = 0.050844906370207;
        
        g(1, 4) = 0.501426509658179;
        g(2, 4) = 0.249286745170910;
        g(3, 4) = 0.249286745170910;
        g(4, 4) = 0.116786275726379;
        
        g(1, 5) = 0.249286745170910;
        g(2, 5) = 0.501426509658179;
        g(3, 5) = 0.249286745170910;
        g(4, 5) = 0.116786275726379;
        
        g(1, 6) = 0.249286745170910;
        g(2, 6) = 0.249286745170910;
        g(3, 6) = 0.501426509658179;
        g(4, 6) = 0.116786275726379;
        
        g(1, 7) = 0.636502499121399;
        g(2, 7) = 0.310352451033785;
        g(3, 7) = 0.053145049844816;
        g(4, 7) = 0.082851075618374;
        
        g(1, 8) = 0.636502499121399;
        g(2, 8) = 0.053145049844816;
        g(3, 8) = 0.310352451033785;
        g(4, 8) = 0.082851075618374;
        
        g(1, 9) = 0.310352451033785;
        g(2, 9) = 0.636502499121399;
        g(3, 9) = 0.053145049844816;
        g(4, 9) = 0.082851075618374;
        
        g(1,10) = 0.053145049844816;
        g(2,10) = 0.636502499121399;
        g(3,10) = 0.310352451033785;
        g(4,10) = 0.082851075618374;
        
        g(1,11) = 0.310352451033785;
        g(2,11) = 0.053145049844816;
        g(3,11) = 0.636502499121399;
        g(4,11) = 0.082851075618374;
        
        g(1,12) = 0.053145049844816;
        g(2,12) = 0.310352451033785;
        g(3,12) = 0.636502499121399;
        g(4,12) = 0.082851075618374;
        
    otherwise
        error(['Gauss integration on triangle: ',num2str(n_g),'-point formula not implemented'])
end

% the first, second and third rows of g are the areal coordinates
% hence, sum(el(1:3,:))=ones
L1=g(1,:);
L2=g(2,:);
L3=g(3,:);
% the fourth row of el is the weight
w=g(4,:);

end
