# subradius_computation

This directory contains all MATLAB scripts' that are necessary to recreate the simulations on the paper ArXiv:??? for the computation of the subradius of a family of real non-negative matrices.

- illustrative_example.m can be run as it is. All three algorithms are applied to the same matrix family (which is also discussed in Section 5.1 of the paper) for different values of M. Every option is customizable such as the target accuracy, maximum number of evaluations, etc.
- real_antinorm.m and matrix_antinorm.m are used to compute, respectively, the polytope antinorm of a vector and of a matrix.
- lsr_comp_standard.m is the standard algorithm, in which the polytope antinorm is fixed in the beginning and does not change at any point.
- adaptive_subradius_comp.m is the adaptive variant in which new vertices are added if and only if they improve the current antinorm.
- adaptive_eigenvectors_subradius_comp.m is another adaptive variant, which further refines the antinorm by adding suitably rescaled eigenvectors. NOTE: as mentioned in the paper, this is the only variant that does not have a theoretical result which guarantees convergence.

The file illustrative_example.m can easily be modified to apply one (or all) algorithm(s) to a different matrix family. In that case, replace A_1 and A_2 accordingly and, since the LSR is likely not known, renormalize the family using, for example, the result of "min(a(A_1),a(A_2))", where a() is the polytope antinorm chosen as the initial one. Note that V_in must be changed to define the new antinorm and, for a family of dxd matrices, it should be a dxp matrix, where the p columns are the vertices.

This project is licensed under the GNU General Public License v3.0. You are free to use, modify, and distribute the code for both non-commercial and commercial purposes. However, the authors provide no warranties and disclaim any liability for the code's use.

Please see the LICENSE file for more details.
