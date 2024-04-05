This is a suite of Matlab (R2012b) subroutines implementing the theory developed in the paper:
Caselli F, Bisegna P. Polar decomposition based corotational framework for triangular shell elements with distributed loads.
International Journal for Numerical Methods in Engineering, 2013
DOI: 10.1002/nme.4528

=========================================================================
As pointed out in Felippa and Haugen, CMAME, 2005, "the key operations of adding and removing rigid body motions can be visualized as a front end filter that lies between the assembler/solver and the element library".

This script reports an example of how the assembler of a general purpose finite element code should use the toolkit.

The toolkit is composed of the following two main routines:

=========================================================================
filter_in, implementing the closed-form formulas derived in Section 2. It is called by the assembler before the core-element routine, and has the purpose of removing rigid body motions. In particular, it transforms the element nodal parameters a into their filtered counterpart bar_a, and the dead [respectively, follower] distributed loads into their pulled-back counterparts, to be passed to the core element. Its dependencies are the routines:

(i) rot2spin, transforming variations of nodal rotation vectors to nodal spins and returning the tensor M (Section 2.3.1);
(ii) hat, computing the rototranslation (hat_R, t) and the related tensors hat_G, hat_P (Section 2.3.2);
(iii) chk, computing the rotation chk_R and the related tensors chk_G, chk_P (Section 2.3.3);
(iv) spin2rot,  transforming filtered nodal spins to variations of filtered nodal rotation vectors and returning the tensor B (Section 2.3.4);

=========================================================================
filter_out, implementing the equations presented in Sections 3 and 4. It is called by the assembler after the core-element routine, and has the purpose of computing the nodal residual vector q and the consistent tangent stiffness tensor K, respectively through the algorithms (97), (98), and (127), (128), by processing the following quantities supplied by the core element:

(i) bar_q_i, bar_K_bar_q_i, nodal forces equivalent to internal forces and material tangent stiffness tensor;
(ii) bar_q_load, bar_Q_load, bar_QM_load, bar_H_load, opposite of nodal forces equivalent to applied loads
(iii) f_load, m_G_load, M_G_load, resultants of applied loads

Dependencies of filter_out are the routines:
(i) D_Bt, returning chk_K_Bt, equation (B.3);
(ii) D_shp_R_Gt, returning hat_K_shp_chk_R_chk_Gt, equation (111), and tld_K_shp_hat_R_hat_Gt, equation (118);
(iii) D_shp_R_Pt, returning hat_K_shp_chk_R_chk_Pt, equation (125), and tld_K_shp_hat_R_hat_Pt, equation (126);
(vi) D_Mt, returning K_Mt, equation (B.1).

========================================================================= 
Computations are performed using coordinates in the element reference frame (e,h,n) which is returned by the routine geometry. Other auxiliary self-explanatory routines are included.

=========================================================================
For the sake of completeness, a small-strain core-element implementation is provided as an example. In particular, a combination of the discrete Kirchhoff plate triangle (DKT) and the optimal ANDES template membrane triangle (OPT) is implemented in the routine DKT_OPT_shell: details are omitted here, since that implementation is standard. Of course, any other small- or finite-strain triangular shell element with the same node and degree-of-freedom configuration can be used in the present corotational framework.

=========================================================================
Authors' e-mail addresses: 
caselli@ing.uniroma2.it (Federica Caselli)
bisegna@uniroma2.it (Paolo Bisegna)

(C) 2010-2013 Paolo Bisegna and Federica Caselli. License: GNU General Public License (GPLv3)
