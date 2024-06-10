# subradius_computation

This directory contains all MATLAB scripts' that are necessary to recreate the simulations on the paper ArXiv:??? for the computation of the subradius of a family of real non-negative matrices.

- illustrative_example.m can be run as it is. All three algorithms are applied to the same matrix family (which is also discussed in Section 5.1 of the paper) for different values of M. Every option is customizable and it is possible to change entirely the matrix family. In that case, some care is required (i.e., remove the renormalization and choose V_in as eye(d), d being the correct dimension).
- real_antinorm.m and matrix_antinorm.m are used to compute, respectively, the polytope antinorm of a vector and of a matrix.
- lsr_comp_standard.m is the standard algorithm, in which the polytope antinorm is fixed in the beginning and does not change at any point.
- adaptive_subradius_comp.m is the adaptive variant in which new vertices are added if and only if they improve the current antinorm.
- adaptive_eigenvectors_subradius_comp.m is another adaptive variant, which further refines the antinorm by adding suitably rescaled eigenvectors. NOTE: as mentioned in the paper, this is the only variant that does not have a theoretical result which guarantees convergence.
