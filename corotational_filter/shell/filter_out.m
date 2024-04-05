% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528

function [ ... % nodal residual vector and consistent tangent stiffness tensor
    q, K] = filter_out( ...
    ... % nodal parameters and filtered nodal parameters
    a, bar_a, ...
    ... % tensor relating delta_tld_a to delta_a
    M, ...
    ... % tensors needed to filter out the rigid rototranslation t,hat_R
    hat_R, T, hat_G, hat_P, ...
    ... % tensors needed to filter out the rigid rotation chk_R
    chk_R, chk_G, chk_P, ...
    ... % tensor relating delta_bar_a to delta_chk_a
    B, ...
    ... % derivative of deformed area with respect to bar_a, divided by deformed area
    b, ...
    ... % nodal forces equivalent to internal forces and material tangent stiffness tensor
    bar_q_i, bar_K_bar_q_i, ...
    ... % opposite of nodal forces equivalent to applied loads
    bar_q_load, bar_Q_load, bar_QM_load, bar_H_load, ...
    ... % resultants of applied loads
    f_load, m_G_load, M_G_load)

% This routine implements the equations presented in Sections 3 and 4. 
% It is called by the assembler after the core-element routine, 
% and has the purpose of computing the nodal residual vector q
% and the consistent tangent stiffness tensor K, 
% respectively through the algorithms (97),(98), and (127),(128), 
% by processing the following quantities supplied by the core element:
%   (i)   bar_q_i, bar_K_bar_q_i, nodal forces equivalent to internal forces and material tangent stiffness tensor
%   (ii)  bar_q_load, bar_Q_load, bar_QM_load, bar_H_load, opposite of nodal forces equivalent to applied loads;
%   (iii) f_load, m_G_load, M_G_load, resultants of applied loads.
% note: bar_QM_load is related to distributed couples, not taken into account in the paper
%
% input:
% 
% a(18x1), bar_a(18x1)      nodal parameters and filtered nodal parameters
% 
% M(18x18)                  tensor relating delta_tld_a to delta_a
%                           delta_tld_a = M*delta_a
% 
% hat_R(3x3), T(3x18), hat_G(3x18), hat_P(18x18) tensors needed to filter out the rigid rototranslation t,hat_R
%                           delta_hat_a = hat_P*shp_hat_R'*delta_tld_a
% 
% chk_R(3x3), chk_G(3x18), chk_P(18x18) tensors needed to filter out the rigid rotation chk_R
%                           delta_chk_a = chk_P*shp_chk_R'*delta_hat_a
% 
% B(18x18)                  tensor relating delta_bar_a to delta_chk_a
%                           delta_bar_a = B*delta_chk_a
% 
% b(1x18)                   derivative of deformed area with respect to bar_a, divided by deformed area
% 
% bar_q_i(18x1), bar_K_bar_q_i(18x18) nodal forces equivalent to internal forces and material tangent stiffness tensor
% 
% bar_q_load(18x2), bar_Q_load(18x3x2), bar_QM_load(18x3x2), bar_H_load(18x18x2) opposite of nodal forces equivalent to applied loads
%                           last index: 1=dead loads; 2=follower loads
% 
% f_load(3x2), m_G_load(3x2), M_G_load(3x3x2) resultants of applied loads
%                           last index: 1=dead loads; 2=follower loads
%
% 
% output:
% 
% q(18x1)                   nodal residual vector (in element reference frame)
% 
% K(18x18)                  consistent tangent stiffness tensor (in element reference frame)
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% virtual work: internal force contribution
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% equation (88)
% q_i*delta_a =
%   int_vol_ref S(E)     : delta_E =  (corotational approach)
%   int_vol_ref S(bar_E) : delta_bar_E =
% bar_q_i*delta_bar_a
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% virtual work: load contribution
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% l_d = dead force per unit reference area
% l_f = follower force per unit deformed area
% m_d = dead couple per unit reference area (not accounted for in the paper)
% m_f = follower couple per unit deformed area (not accounted for in the paper)
% those quantities do not depend on nodal parameters
% l_f and m_f rotate according to R
% overall load:
% equation (89)
% l=l_d+R*l_f*area_def/area_ref
% m=m_d+R*m_f*area_def/area_ref
% equation (90)
% q_l*delta_a =
%   -int_area_ref (l*delta_u+m*delta_omega) =  (corotational approach)
%   -(chk_R*f)*delta_tau -(chk_R*m_G)*delta_hat_theta -(m_G)*delta_chk_theta +bar_q_l*delta_bar_a
% with
%   f=f_d+f_f  (equation (91))
%   m_G=m_G_d+m_G_f  (equation (92))
%   bar_q_l=bar_q_d+bar_q_f  (equation (93))
% and
%   f_d, m_G_d force and resultant moment wrt G of R'*l_d and R'*m_d, i.e.
%   f_d = int_area_ref R'*l_d  (equation (91))
%   m_G_d = int_area_ref [p+bar_u(p)-G]x(R'*l_d)+R'*m_d  (equation (92))
%   bar_q_d = -int_area_ref {(shape functions)'*(R'*l_d; R'*m_d)}  (equation (93))
%   f_f, m_G_f force and resultant moment wrt G of l_f and m_f, rescaled with area_def/area_ref, i.e.
%   f_f = int_area_ref l_f*area_def/area_ref  (equation (91))
%   m_G_f = int_area_ref {[p+bar_u(p)-G]x(l_f)+m_f}*area_def/area_ref  (equation (92))
%   bar_q_f = -int_area_ref (shape functions)'*(l_f; m_f)*area_def/area_ref  (equation (93))
%   shape functions = D_delta_bar_u_D_delta_bar_a
% note that, at the core element level,
% the dead loads are transformed into (R'*l_d; R'*m_d)
% the follower loads are transformed into (l_f; m_f) area_def/area_ref
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% virtual work: adding together contributions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% we compute:
% equation (96)
% q*delta_a = tld_q_l*delta_tld_a + hat_q_l*delta_hat_a + (bar_q_i+bar_q_l)*delta_bar_a
%
% with:
% contribution at delta_bar_a (core element) level:
%   from internal forces
%   bar_q_i  (equation (88))
%   from applied dead and follower loads
%   bar_q_l  (equation (93))
%
% contribution at delta_hat_a level:
% equation (95)
% recalling that delta_chk_theta = chk_G*shp_chk_R'*delta_hat_a:
%   from (-mG*delta_chk_theta)
%   hat_q_l = -shp_chk_R*chk_G'*m_G;
%
% contribution at delta_tld_a level:
% equation (95)
% recalling that delta_tau = T*shp_hat_R'*delta_tld_a, 
% and delta_hat_theta = hat_G*shp_hat_R'*delta_tld_a:
%   from (-chk_R*f*delta_tau)
%   tld_q_l_I = -shp_hat_R*T'*chk_R*f;
%   from (-chk_R*mG*delta_hat_theta)
%   tld_q_l_II = -shp_hat_R*hat_G'*chk_R*m_G;
%   summing:
%   tld_q_l=tld_q_l_I+tld_q_l_II
%
% the virtual work is recast as follows:
% equation (97)
%   tld_q_l*delta_tld_a + hat_q_l*delta_hat_a + (bar_q_i+bar_q_l)*delta_bar_a
% (using delta_bar_a = B*delta_chk_a)
%   tld_q_l*delta_tld_a + hat_q_l*delta_hat_a + chk_q*delta_chk_a =
% (using delta_chk_a = chk_P*shp_chk_R'*delta_hat_a)
%   tld_q_l*delta_tld_a + hat_q*delta_hat_a =
% (using delta_hat_a = hat_P*shp_hat_R'*delta_tld_a)
%   tld_q*delta_tld_a =
% (using delta_tld_a = M*delta_a)
%   q*delta_a
%
% with:
% equation (98)
%   bar_q = bar_q_i + bar_q_l
%   chk_q = B'*bar_q
%   hat_q = hat_q_l + shp_chk_R*chk_P'*chk_q
%   tld_q = tld_q_l + shp_hat_R*hat_P'*hat_q;
%   q = M'*tld_q
%
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% linearization of virtual work
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% equation (102)
% Delta(q*delta_a) =
% Delta((bar_q_i+bar_q_l)*delta_bar_a + hat_q_l*delta_hat_a + tld_q_l*delta_tld_a) =
%   Delta(bar_q_i+bar_q_l)*delta_bar_a + Delta(hat_q_l)*delta_hat_a + Delta(tld_q_l)*delta_tld_a
%   + (bar_q_i+bar_q_l)*Delta(delta_bar_a) + hat_q_l*Delta(delta_hat_a) + tld_q_l*Delta(delta_tld_a)
%
% ... after some algebra ...
%
% equation (127)
%   bar_K[Delta_bar_a]*delta_bar_a + hat_K_l[Delta_hat_a]*delta_hat_a + tld_K_l[Delta_tld_a]*delta_tld_a
%   + bar_q*Delta(delta_bar_a) + hat_q_l*Delta(delta_hat_a) + tld_q_l*Delta(delta_tld_a) =
% (using delta_bar_a = B*delta_chk_a)
%   (B'*bar_K*B+chk_K_Bt(bar_q))[Delta_chk_a]*delta_chk_a + hat_K_l[Delta_hat_a]*delta_hat_a + tld_K_l[Delta_tld_a]*delta_tld_a
%   + chk_q*Delta(delta_chk_a) + hat_q_l*Delta(delta_hat_a) + tld_q_l*Delta(delta_tld_a) =
% (using delta_chk_a = chk_P*shp_chk_R'*delta_hat_a)
%   (shp_chk_R*chk_P'*chk_K*chk_P*shp_chk_R'+hat_K_shp_chk_R_chk_Pt(chk_q)+hat_K_l)[Delta_hat_a]*delta_hat_a + tld_K_l[Delta_tld_a]*delta_tld_a
%   + hat_q*Delta(delta_hat_a) + tld_q_l*Delta(delta_tld_a) =
% (using delta_hat_a = hat_P*shp_hat_R'*delta_tld_a)
%   (shp_hat_R*hat_P'*hat_K*hat_P*shp_hat_R'+tld_K_shp_hat_R_hat_Pt(hat_q)+tld_K_l)[Delta_tld_a]*delta_tld_a
%   + tld_q*Delta(delta_tld_a) =
% (using delta_tld_a = M*delta_a)
%   (M'*tld_K*M+K_Mt(tld_q))[Delta_a]*delta_a =
%   K[Delta_a]*delta_a
%
% with:
% equations (122) and (128)
% bar_K = bar_K_bar_q_i + bar_K_bar_q_l
% chk_K = B'*bar_K*B + chk_K_Bt(bar_q)
% hat_K = shp_chk_R*chk_P'*chk_K*chk_P*shp_chk_R' + hat_K_shp_chk_R_chk_Pt(chk_q) + hat_K_l
% tld_K = shp_hat_R*hat_P'*hat_K*hat_P*shp_hat_R' + tld_K_shp_hat_R_hat_Pt(hat_q) + tld_K_l
% K     = M'*tld_K*M + K_Mt(tld_q)
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

% number of element nodes
n_el_n = 3;

% indices of translational and rotational dofs
[~, ind_rot]=indices(n_el_n);

% nodal rotation vectors
% theta_i(:,i) is the nodal rotation vector of node i
theta_i=a(ind_rot);
% core nodal rotation vectors
% bar_theta_i(:,i) is the core nodal rotation vector of node i
bar_theta_i=bar_a(ind_rot);

% expand rotation matrices
shp_chk_R = shp(chk_R, n_el_n);
shp_hat_R = shp(hat_R, n_el_n);

% useful matrix products
chk_G_shp_chk_Rt=chk_G*shp_chk_R';
chk_P_shp_chk_Rt=chk_P*shp_chk_R';
hat_G_shp_hat_Rt=hat_G*shp_hat_R';
hat_P_shp_hat_Rt=hat_P*shp_hat_R';
% ... and their transposes
shp_chk_R_chk_Gt=chk_G_shp_chk_Rt';
shp_chk_R_chk_Pt=chk_P_shp_chk_Rt';
shp_hat_R_hat_Gt=hat_G_shp_hat_Rt';
shp_hat_R_hat_Pt=hat_P_shp_hat_Rt';

% construct flt_I=blkdiag(I,O,I,O,I,O)
I = eye(3);
O = zeros(3);
flt_I = blkdiag(I,O,I,O,I,O); % n_el_n copies

% extract contribution from dead loads
% opposite of nodal forces equivalent to applied loads
% -int( (displacement shape functions)'*(distributed forces) + (rotation shape functions)'*(distributed couples))
% equation (93)
bar_q_d=bar_q_load(:,1);
% -int( (displacement shape functions)'*spin(distributed forces) )
% equation (105)
bar_Q_d=bar_Q_load(:,:,1);
% -int( (displacement shape functions)'*spin(distributed forces) + (rotation shape functions)'*spin(distributed couples) )
bar_QM_d=bar_QM_load(:,:,1);
% -int( Delta_(displacement shape functions)'*(distributed forces) + Delta_(rotation shape functions)'*(distributed couples) )
bar_H_d=bar_H_load(:,:,1);
% resultant force: int( distributed forces )
% equation (91)
f_d=f_load(:,1);
% resultant moment wrt G: int( spin(position vector wrt G)*(distributed forces) + (distributed couples) )
% equation (92)
m_G_d=m_G_load(:,1);
% int( spin(position vector wrt G)*spin(distributed forces) + spin(distributed couples) )
% equation (113)
M_G_d=M_G_load(:,:,1);

% extract contribution from follower loads
% opposite of nodal forces equivalent to applied loads
% -int( (displacement shape functions)'*(distributed forces) + (rotation shape functions)'*(distributed couples))
% equation (93)
bar_q_f=bar_q_load(:,2);
% -int( (displacement shape functions)'*spin(distributed forces) )
% equation (113)
bar_Q_f=bar_Q_load(:,:,2);
% -int( Delta_(displacement shape functions)'*(distributed forces) + Delta_(rotation shape functions)'*(distributed couples) )
bar_H_f=bar_H_load(:,:,2);
% resultant force: int( distributed forces )
% equation (91)
f_f=f_load(:,2);
% resultant moment wrt G: int( spin(position vector wrt G)*(distributed forces) + (distributed couples) )
% equation (92)
m_G_f=m_G_load(:,2);

% adding together dead and follower load contributions
% equations (93); text between (112) and (113); (91); (92)
bar_q_l = bar_q_d + bar_q_f;
bar_Q_l = bar_Q_d + bar_Q_f;
bar_H_l = bar_H_d + bar_H_f;
f       = f_d     + f_f;
m_G     = m_G_d   + m_G_f;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% virtual work
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% contribution at delta_bar_a (core element) level:
%   from internal forces
%   bar_q_i
%   from applied dead and follower loads
%   bar_q_l
% equation (98)
bar_q = bar_q_i + bar_q_l;

% contribution at delta_chk_a level
%   none
% transformation of bar_q to delta_chk_a level
% equation (98)
chk_q = B'*bar_q;

% contribution at delta_hat_a level,
%   from (-mG*delta_chk_theta)
% equation (95)
hat_q_l = -shp_chk_R_chk_Gt*m_G;
% transformation of chk_q to delta_hat_a level and load contribution at this level
% equation (98)
hat_q = hat_q_l + shp_chk_R_chk_Pt*chk_q;

% contribution at delta_tld_a level:
% equation (95)
tld_q_l = ...
    ... % from (-chk_R*f*delta_tau)
    -shp_hat_R*T'*chk_R*f ...
    ... % from (-chk_R*mG*delta_hat_theta)
    -shp_hat_R_hat_Gt*chk_R*m_G;
% transformation of hat_q to delta_tld_a level and load contribution at this level
% equation (98)
tld_q = tld_q_l + shp_hat_R_hat_Pt*hat_q;

% contribution at delta_a level
%   none
% transformation of tld_q to delta_a level
% element nodal residual vector
% equation (98)
q = M'*tld_q;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% linearization of virtual work
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% contribution from the term involving Delta(bar_q_i)
% equation (103)
%   Delta(bar_q_i)*delta_bar_a = bar_K_bar_q_i[Delta_bar_a]*delta_bar_a
% with:
% bar_K_bar_q_i = D_q_i_D_bar_a

% contribution from the term involving Delta(bar_q_l)
% equation (107)
%   Delta(bar_q_l)*delta_bar_a = bar_K_bar_q_l[Delta_bar_a]*delta_bar_a + hat_K_bar_q_l[Delta_hat_a]*delta_hat_a + tld_K_bar_q_l[Delta_tld_a]*delta_tld_a
% equation (108)
bar_K_bar_q_l = ...
    ... % from (bar_q_f*delta_bar_a)
    bar_q_f*b ...
    ... % from shape function derivative, not taken into account in the paper
    + bar_H_l;
% equation (108): bar_QM_d instead of bar_Q_d takes into account distributed couples, too
hat_K_bar_q_l = ...
    ... % from (bar_q_d*delta_bar_a)
    shp_chk_R_chk_Pt*B'*bar_QM_d*chk_G_shp_chk_Rt;
% equation (108): bar_QM_d instead of bar_Q_d takes into account distributed couples, too
tld_K_bar_q_l = ...
    ... % from (bar_q_d*delta_bar_a)
    shp_hat_R_hat_Pt*shp_chk_R_chk_Pt*B'*bar_QM_d*chk_R'*hat_G_shp_hat_Rt ;

% contribution from the term involving Delta(hat_q_l)
% equation (114)
%   Delta(hat_q_l)*delta_hat_a = hat_K_hat_q_l[Delta_hat_a]*delta_hat_a + tld_K_hat_q_l[Delta_tld_a]*delta_tld_a
% equation (115)
hat_K_hat_q_l = ...
    ... % from (-mG_d*delta_chk_theta) and (-mG_f*delta_chk_theta)
    - D_shp_R_Gt(m_G, shp_chk_R, chk_G, chk_P, flt_I, n_el_n) ...
    + shp_chk_R_chk_Gt* ( ...
    + (bar_Q_l'- m_G_f*b)*B*chk_P_shp_chk_Rt ...
    - M_G_d*chk_G_shp_chk_Rt ...
    );
% equation (115)
tld_K_hat_q_l = ...
    ...%   from (-mG_d*delta_chk_theta)
    - shp_hat_R_hat_Pt*shp_chk_R_chk_Gt*M_G_d*chk_R'*hat_G_shp_hat_Rt;


% contribution from the term involving Delta(tld_q_l)
% equation (119)
%   Delta(tld_q_l)*delta_tld_a = tld_K_tld_q_l[Delta_tld_a]*delta_tld_a
% equation (120)
tld_K_tld_q_l = ...
    ... % from (-chk_R*f_f*delta_tau)
    T'*hat_R*chk_R*( ...
    spin(f_f)*(chk_R'*hat_G+chk_G_shp_chk_Rt*hat_P) ...
    -(f_f*b)*B*chk_P_shp_chk_Rt*hat_P ...
    )*shp_hat_R' ...
    ... % from (-chk_R*mG_d*delta_hat_theta) and (-chk_R*mG_f*delta_hat_theta)
    - D_shp_R_Gt(chk_R*m_G, shp_hat_R, hat_G, hat_P, flt_I, n_el_n) ...
    + shp_hat_R_hat_Gt*chk_R*( ...
    + spin(m_G)*chk_G_shp_chk_Rt*hat_P  ...
    + (bar_Q_l'- m_G_f*b)*B*chk_P_shp_chk_Rt*hat_P ...
    - M_G_d*(chk_R'*hat_G+chk_G_shp_chk_Rt*hat_P) ...
    )*shp_hat_R';


% adding together load contribution
% equation (121)
%   Delta(bar_q_i+bar_q_l)*delta_bar_a + Delta(hat_q_l)*delta_hat_a + Delta(tld_q_l)*delta_tld_a =
%   (bar_K_i+bar_K_l)[Delta_bar_a]*delta_bar_a + hat_K_l[Delta_hat_a]*delta_hat_a + tld_K_l[Delta_tld_a]*delta_tld_a
% with:
% equation (122)
% internal force contribution
bar_K_i=bar_K_bar_q_i;
% load contribution
bar_K_l = bar_K_bar_q_l;
hat_K_l = hat_K_bar_q_l + hat_K_hat_q_l;
tld_K_l = tld_K_bar_q_l + tld_K_hat_q_l + tld_K_tld_q_l;


% contributions from the terms (bar_q_i+bar_q_l)*Delta(delta_bar_a) + hat_q_l*Delta(delta_hat_a) + tld_q_l*Delta(delta_tld_a)
% equations (102), (123)

% contribution at [Delta_chk_a]*delta_chk_a level
% equations (124), (B.3)
%   from bar_q*Delta(delta_bar_a), contribution from B
chk_K_Bt = D_Bt(bar_q, bar_theta_i, n_el_n);

% contribution at [Delta_hat_a]*delta_hat_a level
% equations (124), (125)
%   from chk_q*Delta(delta_chk_a), contribution from chk_P*shp_chk_R'
hat_K_shp_chk_R_chk_Pt = D_shp_R_Pt(chk_q, shp_chk_R, chk_G, chk_P, flt_I, n_el_n);

% contribution at [Delta_tld_a]*delta_tld_a level
% equations (124), (126)
%   from hat_q*Delta(delta_hat_a), contribution from shp_hat_R*hat_P'
tld_K_shp_hat_R_hat_Pt=D_shp_R_Pt(hat_q, shp_hat_R, hat_G, hat_P, flt_I, n_el_n);

% contribution at [Delta_a]*delta_a level
% equations (124), (B.1)
%   from tld_q*Delta(delta_tld_a), contribution from M
K_Mt=D_Mt(tld_q, theta_i, n_el_n);


% computation of the consistent tangent stiffness tensor
% equation (127)

% [Delta_bar_a]*delta_bar_a (core element) level:
% equation (122)
bar_K = bar_K_i + bar_K_l;

% [Delta_chk_a]*delta_chk_a level
% equation (128)
chk_K = ...
    ... % transformation of bar_K to [Delta_chk_a]*delta_chk_a level
    + B'*bar_K*B ...
    ... % contribution from bar_q*Delta(delta_bar_a), contribution from B
    + chk_K_Bt;


% [Delta_hat_a]*delta_hat_a level
% equation (128)
hat_K = ...
    ... % transformation of chk_K to [Delta_hat_a]*delta_hat_a level
    + shp_chk_R_chk_Pt*chk_K*chk_P_shp_chk_Rt ...
    ... % contribution from chk_q*Delta(delta_chk_a), contribution from chk_P*shp_chk_R'
    + hat_K_shp_chk_R_chk_Pt ...
    ... % load contribution
    + hat_K_l;


% [Delta_tld_a]*delta_tld_a level
% equation (128)
tld_K = ...
    ... % transformation of hat_K to [Delta_tld_a]*delta_tld_a level
    + shp_hat_R_hat_Pt*hat_K*hat_P_shp_hat_Rt ...
    ... % contribution from hat_q*Delta(delta_hat_a), contribution from shp_hat_R*hat_P'
    + tld_K_shp_hat_R_hat_Pt ...
    ... % load contribution
    + tld_K_l;


% [Delta_a]*delta_a level
% element consistent tangent stiffness tensor
% equation (128)
K = ...
    ... % transformation of tld_K to [Delta_a]*delta_a level
    + M'*tld_K*M ...
    ... % contribution from tld_q*Delta(delta_tld_a), contribution from M
    + K_Mt;

end
