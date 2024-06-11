# Subradius Computation

This project contains MATLAB scripts for computing the subradius of a family of real non-negative matrices, as described in the paper **ArXiv:???**. The algorithms provided here can be used to recreate all simulations presented in both the main paper and the supplementary material (or appendix, in the ArXiv version).

## Contents

- `illustrative_example.m`: A customizable script that applies all three algorithms to the same matrix family discussed in Section 5.1 of the paper, allowing for different values of M, target accuracy, maximum number of evaluations, etc.
- `random_matrices_simulation.m`: A customizable script that applies the adaptive algorithms, following the strategy discussed in Section 6.2 of the paper, to randomly generated families of either full or sparse matrices. Users can easily customize options such as tested dimensions, sparsity density, maximum number of iterations, initial antinorm, and the specific adaptive algorithm used. This script can be used to make comparisons with the data in Section 7 as they were obtained through a similar strategy.
- `real_antinorm.m`: Computes the polytope antinorm of a vector.
- `matrix_antinorm.m`: Computes the polytope antinorm of a matrix.
- `lsr_comp_standard.m`: The standard algorithm, in which the polytope antinorm is fixed at the beginning and does not change.
- `adaptive_subradius_comp.m`: An adaptive variant that adds new vertices only if they improve the current antinorm.
- `adaptive_eigenvectors_subradius_comp.m`: Another adaptive variant that further refines the antinorm by adding suitably rescaled eigenvectors. Note: As mentioned in the paper, this is the only variant without a theoretical result guaranteeing convergence.

## Usage

Both scripts `illustrative_example.m` and `random_matrices_simulation.m` can be run as they are and, if necessary, the customizable options can be easily modified following the comments in those files.

To apply the algorithms to a different matrix family, modify the `illustrative_example.m` file:

1. Replace `A_1` and `A_2` with the desired matrices.
2. Since the LSR is likely not known, renormalize the family using, for example, the result of `min(a(A_1),a(A_2))`, where `a()` is the chosen initial polytope antinorm. Alternatively, rescale the family using a preliminary bound as in `random_matrices_simulation.m`.
3. Modify `V_in` to define the new antinorm. For a family of dxd matrices, `V_in` should be a dxp matrix, where the p columns are the vertices. For example, `V_in = eye(d)` is the 1-antinorm in the dxd setting.

To apply the perturbation theory to a matrix family of dxd matrices:

1. Introduce the perturbation matrix `Delta = rand(d,d * m)`, were m is the number of matrices in the family (e.g., m = 2 above). Next, rescale the matrix by setting `Delta = Delta/norm(Delta,'fro')`.
2. Let `eps` be the size of the perturbation (e.g., between 10e-3 and 10e-7) and define the perturbed family `A_p = A + eps * Delta`.
3. Apply the adaptive algorithms following either `illustrative_example.m` or `random_matrices_simulation.m` depending on the context to the perturbed family A_p. Possibly, do so for several values of `eps` for a more complete picture.
4. **Optional step:** Verify if the s.l.p. found for the perturbed family when `eps` is small is also a s.l.p. for the initial family when `eps = 0` using, for example, the polytopic algorithm.

## License

This project is licensed under the GNU General Public License v3.0. See the LICENSE file for more details.
