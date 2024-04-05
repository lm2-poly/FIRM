% This is a suite of Matlab (R2012b) subroutines
% implementing the theory developed in the paper:
% Caselli F, Bisegna P. Polar decomposition based corotational framework
% for triangular shell elements with distributed loads.
% International Journal for Numerical Methods in Engineering, 2013
% DOI: 10.1002/nme.4528
% Equation and Section numbers below refer to this paper. 
% 
% As pointed out in Felippa and Haugen, CMAME, 2005, "the key operations of
% adding and removing rigid body motions can be visualized as a front end
% filter that lies between the assembler/solver and the element library".
%
% This script reports an example of how the assembler of a general purpose
% finite element code should use the toolkit.
%
% The toolkit is composed of the following two main routines:
%
% filter_in, implementing the closed-form formulas
%     derived in Section 2. It is called by the assembler before the
%     core-element routine, and has the purpose of removing rigid body
%     motions. In particular, it transforms the element nodal parameters a
%     into their filtered counterpart bar_a, and the dead [respectively,
%     follower] distributed loads into their pulled-back counterparts, to
%     be passed to the core element. Its dependencies are the routines:
%
%     (i) rot2spin, transforming variations of nodal rotation vectors to
%     nodal spins and returning the tensor M (Section 2.3.1);
%
%     (ii) hat, computing the rototranslation (hat_R, t) and the related
%     tensors hat_G, hat_P (Section 2.3.2);
%
%     (iii) chk, computing the rotation chk_R and the related tensors
%     chk_G, chk_P (Section 2.3.3);
%
%     (iv) spin2rot,  transforming filtered nodal spins to variations of
%     filtered nodal rotation vectors and returning the tensor B (Section
%     2.3.4);
%
% filter_out, implementing the equations presented in
%     Sections 3 and 4. It is called by the assembler after the
%     core-element routine, and has the purpose of computing the nodal
%     residual vector q and the consistent tangent stiffness tensor K,
%     respectively through the algorithms (97), (98), and (127), (128), by
%     processing the following quantities supplied by the core element:
%
%     (i) bar_q_i, bar_K_bar_q_i, nodal forces equivalent to internal
%     forces and material tangent stiffness tensor;
%
%     (ii) bar_q_load, bar_Q_load, bar_QM_load, bar_H_load, opposite of
%     nodal forces equivalent to applied loads
%
%     (iii) f_load, m_G_load, M_G_load, resultants of applied loads
%
%     Dependencies of filter_out are the routines:
%
%     (i) D_Bt, returning chk_K_Bt, equation (B.3);
%
%     (ii) D_shp_R_Gt, returning hat_K_shp_chk_R_chk_Gt, equation (111),
%     and tld_K_shp_hat_R_hat_Gt, equation (118);
%
%     (iii) D_shp_R_Pt, returning hat_K_shp_chk_R_chk_Pt, equation (125),
%     and tld_K_shp_hat_R_hat_Pt, equation (126);
%
%     (vi) D_Mt, returning K_Mt, equation (B.1).
%
%
% Computations are performed using coordinates in the element reference
% frame (e,h,n) which is returned by the routine geometry. Other auxiliary
% self-explanatory routines are included.
%
% For the sake of completeness, a small-strain core-element implementation
% is provided as an example. In particular, a combination of the discrete
% Kirchhoff plate triangle (DKT) and the optimal ANDES template membrane
% triangle (OPT) is implemented in the routine DKT_OPT_shell: details are
% omitted here, since that implementation is standard. Of course, any other
% small- or finite-strain triangular shell element with the same node and
% degree-of-freedom configuration can be used in the present corotational
% framework.
%
% Authors' e-mail addresses: 
% caselli@ing.uniroma2.it (Federica Caselli)
% bisegna@uniroma2.it (Paolo Bisegna)
% 
% (C) 2010-2013 Paolo Bisegna and Federica Caselli. License: GNU General Public License (GPLv3)

% clear
clear all
close all
clc
format long e


% The following instructions are supposed to be part of the assembler
% routine, inside a cycle on shell elements

% nodal coordinates of reference triangle (in global frame)
% coor_el(3x3)
% coor_el(:,i) contains the coordinates of node V_i
coor_el = [ ...
    -6.3e+00   +8.2e-01  -7.5e+00; ...
    +1.5e+01   +1.6e+01  +2.4e+01; ...
    +2.8e+01   +1.7e+01  +2.3e+01];

% current nodal parameters (in global frame)
% a_el(18x1)
% a_el=(u_1;theta_1;u_2;theta_2;u_3;theta_3)
a_el = [ ...
    ... % displacement of node 1
    1.7e+001; ...
    -8.9e+000; ...
    5.6e+000; ...
    ... % rotation vector of node 1
    8.4e-001; ...
    -5.2e-001; ...
    -4.0e-001; ...
    ... % displacement of node 2
    2.1e+001; ...
    -4.6e+000; ...
    1.5e+001; ...
    ... % rotation vector of node 2
    8.5e-001; ...
    -4.6e-001; ...
    -4.1e-001; ...
    ... % displacement of node 3
    2.4e+001; ...
    -7.5e+000; ...
    1.5e+001; ...
    ... % rotation vector of node 3
    8.4e-001; ...
    -5.2e-001; ...
    -3.3e-001];

% dead distributed loads per unit reference area
% d_load_el(6x3)
% first index -> load component (qx,qy,qz,mx,my,mz) in GLOBAL reference frame
% second index -> node (1, 2, 3): linear variation over the triangle is implied:
% for a uniform load the three columns of d_load_el are equal
d_load_el = [ ...
    5.0e-01   2.8e-01   8.6e-01; ...
    2.2e-01   1.4e-02   3.2e-01; ...
    2.1e-01   8.1e-01   1.3e-01; ...
    6.8e-01   7.2e-01   7.8e-01; ...
    6.8e-01   1.0e-01   2.2e-02; ...
    3.8e-01   8.5e-02   9.8e-01];

% follower distributed loads per unit deformed area
% f_load_el(6x3)
% first index -> load component (qx,qy,qz,mx,my,mz) in ELEMENT reference frame (e,h,n)
% second index -> node (1, 2, 3): linear variation over the triangle is implied:
% for a uniform load the three columns of f_load_el are equal
% follower loads rotate according to the rotation tensor R and are rescaled by area_def/area
f_load_el = [ ...
    9.1e-01   1.6e-02   6.4e-01; ...
    3.1e-01   8.0e-01   4.9e-02; ...
    6.0e-01   8.9e-01   2.6e-02; ...
    6.0e-01   2.1e-01   1.7e-01; ...
    9.5e-01   6.7e-01   6.1e-01; ...
    6.3e-02   4.8e-01   8.6e-01];

% type of constitutive law
type='linear_elastic_isotropic_plane_stress';
% Young modulus
E=1e2;
% Poisson's ratio
nu=0.3;
% constitutive parameters
param=struct('E', E, 'nu',nu);
% constitutive properties
constitutive=struct('type',type, 'param',param);

% thickness
th=0.5;

% n_g: number of Gauss points used for numeric integration over triangle
%      admissible values are:
%       1:       1-point gauss integration
%       3:       3-point gauss integration
%       4:       4-point gauss integration
%       7:       7-point gauss integration
%      12:      12-point gauss integration
n_g=7;

% material
mat_el=struct('n_g',n_g, 'constitutive',constitutive, 'th',th);

% Call geometry:
% This routine computes geometrical features of reference and deformed triangle
% (more precisely, the latter is the triangle joining corrent nodal positions);
% moreover, it prepares computations in element reference frame
[ ...
    ... % reference triangle: centroid G and element reference frame (e,h,n) (in global frame)
    ~, ref2glob, ...
    ... % reference triangle: nodal coordinates (in element reference frame), area, side lengths
    coor, area, l23, l31, l12, ...
    ... % deformed triangle: centroid G_def and attatched frame (e_def,h_def,n_def) (in global frame)
    ~, ~, ...
    ... % deformed triangle: nodal coordinates (in frame attatched to deformed triangle), area, side lengths
    ~, area_def, l23_def, l31_def, l12_def, ...
    ... % deformed triangle: attatched frame (in element reference frame)
    e_def, h_def, n_def, ...
    ... % nodal parameters (in element reference frame)
    a, ...
    ... % dead and follower loads (in element reference frame)
    d_load, f_load]=geometry( ...
    ... % nodal coordinates and nodal parameters (in global frame)
    coor_el, a_el, ...
    ... % dead and follower loads (the former in global frame; the latter in element reference frame)
    d_load_el, f_load_el);

% Call filter_in:
% This routine implements the closed-form formulas derived in Section 2.
% It is called by the assembler before the core-element routine,
% and has the purpose of removing rigid body motions.
% In particular, it transforms the element nodal parameters a
% into their filtered counterpart bar_a,
% and the dead [respectively, follower] distributed loads into their pulled-back counterparts
% R'l_d [respectively, l_f*area_def/area, to be passed to the core element.
[ ... % filtered nodal parameters (to be passed to the core element)
    bar_a, ...
    ... % tensor relating delta_tld_a to delta_a
    M, ...
    ... % tensors needed to filter out the rigid rototranslation t,hat_R
    ~, hat_R, T, hat_G, hat_P, ...
    ... % tensors needed to filter out the rigid rotation chk_R
    chk_R, chk_G, chk_P, ...
    ... % tensor relating delta_bar_a to delta_chk_a
    B, ...
    ... % derivative of deformed area with respect to bar_a, divided by deformed area
    b, ...
    ... % number of load conditions and loads
    n_load, load]=filter_in( ...
    ... % nodal coordinates of reference triangle
    coor, ...
    ... % nodal parameters
    a, ...
    ... % frame attached to deformed triangle
    e_def, h_def, n_def, ...
    ... % area and side lengths of reference triangle
    area, l23, l31, l12, ...
    ... % area and side lengths of deformed triangle
    area_def, l23_def, l31_def, l12_def, ...
    ... % dead and follower loads
    d_load, f_load);

% Call the core-element routine:
% This function implements the triangular shell element: OPT/CST membrane + DKT plate
% and returns:
%   (i)   bar_q_i, bar_K_bar_q_i, nodal forces equivalent to internal forces and material tangent stiffness tensor
%   (ii)  bar_q_load, bar_Q_load, bar_QM_load, bar_H_load, opposite of nodal forces equivalent to applied loads;
%   (iii) f_load, m_G_load, M_G_load, resultants of applied loads.
[ ... % material tangent stiffness matrix
    bar_K_bar_q_i, ...
    ... % nodal forces equivalent to internal forces
    bar_q_i, ...
    ... % opposite of nodal forces equivalent to applied loads
    bar_q_load, bar_Q_load, bar_QM_load, bar_H_load, ...
    ... % resultants of applied loads
    f_load, m_G_load, M_G_load]= DKT_OPT_shell( ...
    ... % nodal coordinates and nodal parameters
    coor, bar_a, ...
    ... % number of load conditions and applied loads
    n_load, load, ...
    ... % material parameters
    mat_el);


% Call filter_out:
% This routine implements the equations presented in Sections 3 and 4.
% It is called by the assembler after the core-element routine,
% and has the purpose of computing the nodal residual vector q
% and the consistent tangent stiffness tensor K,
% respectively through the algorithms (97),(98), and (127),(128),
% by processing the quantities supplied by the core element
[ ... % nodal residual vector and consistent tangent stiffness tensor
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
    f_load, m_G_load, M_G_load);


% The nodal residual vector and consistent tangent stiffness tensor
% supplied by filter_out are given in element reference frame. Before
% being assembled, they must be converted to global frame:
shp_ref2glob=shp(ref2glob,3);
q=shp_ref2glob*q;
K=shp_ref2glob*K*shp_ref2glob';

% display results
disp('This is a suite of Matlab (R2012b) subroutines')
disp('implementing the theory developed in the paper:')
disp(' ')
disp('Caselli F, Bisegna P. Polar decomposition based corotational framework')
disp('for triangular shell elements with distributed loads.')
disp('International Journal for Numerical Methods in Engineering, 2013')
disp('DOI: 10.1002/nme.4528')
disp(' ')
disp('(C) License: GNU General Public License (GPLv3)')
disp(' ')
disp('type ''help example.m'' for documentation')
disp(' ')
disp('the following results should coincide with those pasted at the end of the file example.m')
disp(' ')

% display the nodal residual vector
disp('nodal residual vector')
disp(q)

% display the consistent tangent stiffness tensor
disp('consistent tangent stiffness tensor')
disp(K)

% With the data given above, the following output should be obtained:
%
% nodal residual vector
%     -5.658755322613256e+01
%     -5.412735920083155e+01
%     -3.691534006516709e+01
%     -1.266986427434942e+02
%      3.293625690285118e+01
%     -9.147197172666492e+01
%     -5.826979190666135e+01
%      1.387332708520952e+01
%      2.135698104807110e+01
%     -3.876351853795924e+01
%     -4.256887436633357e+01
%     -2.516435898340438e+01
%      4.400021793108556e+01
%      1.987183653690759e+01
%     -3.333003069425789e+01
%      2.001105640630530e+02
%     -1.003330297866432e+01
%      1.878617277898143e+02
%
% consistent tangent stiffness tensor
%   Columns 1 through 4
%
%      2.585461509101828e+01     1.257592124322431e+01    -1.399908960333528e+00    -1.298614928421175e+01
%      1.124843168374190e+01     2.554783806251054e+01     7.172849968648572e+00     4.468844256027753e+01
%     -7.256841498461659e-01     8.249465354071663e+00     1.164089700232583e+01     2.898756088579168e+01
%     -1.230517184675056e+01     4.650611113709354e+01     3.019885477525300e+01     2.411705601739619e+02
%     -1.485382931700091e+00    -3.411978189692419e+00    -5.226988971995962e+00    -3.275410230360087e+01
%     -1.339919761759821e+01     4.972717971924649e+01     4.298521720623889e+01     2.875539347363401e+02
%     -1.272580614893585e+01    -5.898611161861744e+00     2.476286329368441e+00     1.797512359443302e+01
%     -4.604295666520615e+00    -1.553614782585589e+00     2.056980989177726e-02    -7.685653300062167e-01
%      3.965982444268106e+00     1.225681451168564e+00     3.462811378670330e-01     6.273291148217255e+00
%      1.658890159141429e+01     1.845963164961184e+01     1.384780405362335e+00     1.122631320799967e+00
%     -1.835256723598187e+00     1.911568249038777e+00    -2.973564965496113e+00     1.063333118340824e-01
%      1.979259008590063e+01     2.161994112298165e+01     3.999206238340211e+00     8.692643842659065e-01
%     -1.032048382121621e+01    -4.720395915869998e+00    -1.932136184286382e+00    -4.988974310221283e+00
%     -7.896284756171813e+00    -2.132962095533161e+01    -8.950514403188638e+00    -4.391987723027130e+01
%     -1.640241167027540e+00    -8.627198685681606e+00    -9.388174387551668e+00    -3.526085203400893e+01
%     -3.920730092173645e+01    -3.748300264892658e+01     1.907194287902758e+00    -1.012747136612102e+01
%      8.504666760678266e+00    -2.645616372056387e+00     7.417320046378669e+00     1.039766787327020e+00
%     -5.183593767728608e+01    -3.621593380194435e+01    -1.085481881076203e+01    -1.326957455881765e+01
%
%   Columns 5 through 8
%
%      1.037250111461958e-01    -1.091722795738996e+01    -1.455841221146519e+01    -4.165086795115109e+00
%     -3.402274469096977e+00     4.771347547316768e+01    -6.618918111428412e+00    -2.790326211820105e+00
%     -3.755913389046590e+00     4.434199770673402e+01     1.437943963817369e+00     1.204132269014549e+00
%     -3.275410230360275e+01     2.875539347363471e+02     1.447768534734122e+01     2.746220896364333e-01
%     -6.250149524515303e+00    -3.299436377843295e+01    -1.680915328868589e+00     3.588331823421922e-01
%     -3.299436377843205e+01     3.302902312420623e+02     1.549009149645841e+01     8.092200796198540e+00
%     -3.584029804212706e+00     1.589213340922158e+01     1.691461676002080e+01     2.498393610695580e+00
%      1.361354304498351e+00     8.695788102636476e+00     1.861200739544628e+00     6.283613409780866e+00
%     -3.878454168615646e+00    -5.421329029810119e+00    -1.222402504116123e+01     8.223938771532815e+00
%     -1.817792746764589e-01     8.508033266391279e-01    -1.487286012379449e+01    -4.380080364018888e+01
%      5.014491589833167e-01    -2.291873254815189e-01    -7.446086817743272e-01     4.716695684097656e+00
%     -6.761057223385603e-03     1.444597189024299e+00    -1.894822988355438e+01    -5.286237082811603e+01
%      3.480304793066515e+00    -4.974905451831626e+00    -4.754186243398364e+00     2.280948064120209e+00
%      2.040920164598624e+00    -5.640926357580413e+01     3.675108987452973e+00    -4.029416586970224e+00
%      7.634367557662237e+00    -3.892066867692390e+01     6.999111298475857e+00    -8.001181784764466e+00
%      1.074079398966138e+00    -1.312152541911356e+01    -1.793492030028807e+01     6.365116121878726e+00
%      2.363615370326231e-01     1.394310521735432e+00     4.001650351601076e+00    -1.998055584699525e+00
%      1.444728331894242e+00    -1.559183284227138e+01    -1.490033279140369e+01    -4.926999923364515e+00
%
%   Columns 9 through 12
%
%      4.904580580407676e+00     1.700538534086008e+01    -4.325212100692037e+00     2.048361839541826e+01
%      1.063756028802989e+00     2.064007808986808e+01     1.431133045401053e+00     2.450489030827380e+01
%     -1.041779846703935e+00     4.653648869207022e-01    -4.882643400563670e+00     4.081343475864271e+00
%      5.287194679218223e+00     1.122631320799976e+00     1.063333118340812e-01     8.692643842659178e-01
%     -3.915243117732668e+00    -1.817792746764597e-01     5.014491589833168e-01    -6.761057223386446e-03
%     -7.181524115814786e+00     8.508033266391391e-01    -2.291873254815201e-01     1.444597189024312e+00
%     -9.555012825720628e+00    -1.692951828249920e+01     2.035954522658013e+00    -2.271274023561201e+01
%      7.943130004554316e+00    -4.337115618210785e+01     3.779953488102301e+00    -5.201000483206231e+01
%      3.959250309781371e+00    -1.879125124411092e+01     2.634213010000941e+00    -2.669167838523094e+01
%     -1.966549025704692e+01     1.767620934801562e+02    -1.615672308664612e+01     2.077259919870903e+02
%      3.925851484018850e+00    -1.615672308664132e+01    -2.882119230962807e-01    -2.443065354520342e+01
%     -2.684445854921274e+01     2.077259919871002e+02    -2.443065354519651e+01     2.514220034199336e+02
%      8.806479948251511e+00    -7.586705836087582e-02     2.289257578034025e+00     2.229121840193765e+00
%     -8.665916872641688e+00     2.273107809223978e+01    -5.211086533503346e+00     2.750511452378850e+01
%     -4.063102493126808e+00     1.832588635719022e+01     2.248430390562731e+00     2.261033490936667e+01
%     -1.675708076723568e+00     4.194858763465994e+00    -5.423666017143720e-01     4.634827801525464e+00
%      3.237741805484726e+00    -2.803184583713099e-01     5.773269809541541e-01    -5.651308130081483e-01
%      2.095456206631221e+01     4.723031184738557e+00    -3.726812311825296e-01     6.105381517965542e+00
%
%   Columns 13 through 16
%
%     -1.129620287955309e+01    -8.410834448109203e+00    -3.504671620074153e+00    -3.912441286421261e+01
%     -4.629513572313490e+00    -2.275751185069044e+01    -8.236605997451566e+00    -4.036731886865336e+01
%     -7.122598139712055e-01    -9.453597623086216e+00    -1.059911715562190e+01     5.234665651604592e+00
%     -2.172513500590664e+00    -4.678073322672998e+01    -3.548604945447124e+01    -1.012747136612104e+01
%      3.166298260568667e+00     3.053145007350242e+00     9.142232089728633e+00     1.074079398966138e+00
%     -2.090893878860214e+00    -5.781938051544504e+01    -3.580369309042411e+01    -1.312152541911358e+01
%     -4.188810611084951e+00     3.400217551166164e+00     7.078726496352192e+00    -1.753978248742767e+01
%      2.743094926975988e+00    -4.729998627195279e+00    -7.963699814446091e+00     7.074054049148096e+00
%      8.258042596893128e+00    -9.449620222701377e+00    -4.305531447648407e+00    -5.039130393596467e+00
%     -1.716041467619805e+00     2.534117199057702e+01     1.828070985168458e+01     4.194858763465993e+00
%      2.579865405372521e+00    -6.628263933136417e+00    -9.522865185227409e-01    -5.423666017143720e-01
%     -8.443602023462660e-01     3.124242970513437e+01     2.284525231087252e+01     4.634827801525464e+00
%      1.507467006461458e+01     2.439447851749788e+00    -6.874343763965128e+00     5.666419535164029e+01
%      4.221175768718840e+00     2.535903754230185e+01     1.761643127583033e+01     3.329326481950527e+01
%     -5.358870131448319e+00     1.662838047044608e+01     1.345127688067848e+01    -1.955352580081256e-01
%      5.714222122202452e+01     3.111788652704784e+01    -2.314862111791989e-01     2.205336657994970e+02
%     -1.250631711227934e+01     4.643671956755911e+00    -1.065506185186338e+01    -1.541910921699153e+01
%      6.673627046868977e+01     4.114293372530885e+01    -1.009974325555019e+01     2.829202352429043e+02
%
%   Columns 17 through 18
%
%      8.922058206570810e+00    -5.509616913117124e+01
%     -2.462163793326862e+00    -3.955737283753214e+01
%      7.989498176605718e+00    -9.622651018454643e+00
%      1.039766787327023e+00    -1.326957455881767e+01
%      2.363615370326219e-01     1.444728331894243e+00
%      1.394310521735433e+00    -1.559183284227141e+01
%      3.076203677758481e+00    -1.091894395080532e+01
%     -1.561553396354197e+00    -5.979933357740008e+00
%      2.905199289445307e+00     1.632956922785538e+01
%     -2.803184583713105e-01     4.723031184738555e+00
%      5.773269809541546e-01    -3.726812311825289e-01
%     -5.651308130081483e-01     6.105381517965540e+00
%     -1.199826188432930e+01     6.601511308197655e+01
%      4.023717189681053e+00     4.553730619527215e+01
%     -1.089469746605103e+01    -6.706918209400747e+00
%     -1.541910921698758e+01     2.829202352428746e+02
%      1.606549327184124e+01    -2.228107479317421e+01
%     -2.228107479317453e+01     3.868424154092350e+02
