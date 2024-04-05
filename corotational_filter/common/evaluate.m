% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function f=evaluate(funct,x)

% This routine evaluates useful functions having a removable singularity at x=0
% 
% input:
%
% x                         real number, near 0
% 
% funct                     string, specifying the function to evaluate
% 
% presently implemented:    'sin(x)/x'
%                           '(1-cos(x))/x^2'
%                           '(x-sin(x))/x^3'
%                           '(x*sin(x)+2*cos(x)-2)/x^4'
%                           '(3*sin(x)-x*cos(x)-2*x)/x^5'
%                           'asin(x)/x'
%                           'eta'
%                           'mu'
% 
% output:
% 
% f                         function value
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

if abs(x)>1e-2
    % if abs(x) is large enough, use Matlab built-in functions
    switch funct
        case 'sin(x)/x'
            f=sin(x)/x;
        case '(1-cos(x))/x^2'
            f=(1-cos(x))/x^2;
        case '(x-sin(x))/x^3'
            f=(x-sin(x))/x^3;
        case '(x*sin(x)+2*cos(x)-2)/x^4'
            f=(x*sin(x)+2*cos(x)-2)/x^4;
        case '(3*sin(x)-x*cos(x)-2*x)/x^5'
            f=(3*sin(x)-x*cos(x)-2*x)/x^5;
        case 'asin(x)/x'
            f=asin(x)/x;
        case 'eta'
            x2=x/2;
            f=(1-x2*cot(x2))/x^2;
        case 'mu'
            f=(x^2+4*cos(x)+x*sin(x)-4)/(4*x^4*(sin(x/2))^2);
        otherwise
            error(['unknown function ',funct])
    end
else
    % else use Taylor expansion around x=0, up to O(x^8)
    xsq=x*x;
    switch funct
        case 'sin(x)/x'
            % 1-(1/6)*x^2+(1/120)*x^4-(1/5040)*x^6
            f=1-xsq/6*(1-xsq/20*(1-xsq/42));
        case '(1-cos(x))/x^2'
            % 1/2-(1/24)*x^2+(1/720)*x^4-(1/40320)*x^6
            f=(1-xsq/12*(1-xsq/30*(1-xsq/56)))/2;
        case '(x-sin(x))/x^3'
            % 1/6-(1/120)*x^2+(1/5040)*x^4-(1/362880)*x^6
            f=(1-xsq/20*(1-xsq/42*(1-xsq/72)))/6;
        case '(x*sin(x)+2*cos(x)-2)/x^4'
            % -1/12+(1/180)*x^2-(1/6720)*x^4+(1/453600)*x^6
            f=-(1-xsq/15*(1-xsq*3/112*(1-xsq*2/135)))/12;
        case '(3*sin(x)-x*cos(x)-2*x)/x^5'
            % -1/60+(1/1260)*x^2-(1/60480)*x^4+(1/4989600)*x^6
            f=-(1-xsq/21*(1-xsq/48*(1-xsq*2/165)))/60;
        case 'asin(x)/x'
            % 1+(1/6)*x^2+(3/40)*x^4+(5/112)*x^6
            f=(1+xsq/6*(1+xsq*9/20*(1+25/42*xsq)));
        case 'eta'
            % 1/12+(1/720)*x^2+(1/30240)*x^4+(1/1209600)*x^6
            f=(1+xsq/60*(1+xsq/42*(1+xsq/40)))/12;
        case 'mu'
            % 1/360+(1/7560)*x^2+(1/201600)*x^4+(1/5987520)*x^6
            f=(1+xsq/21*(1+xsq*3/80*(1+xsq*10/297)))/360;
        otherwise
            error(['unknown function ',funct])
    end
end
end
