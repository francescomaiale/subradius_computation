# Subradius Computation

This project contains MATLAB scripts for computing the subradius of a family of real non-negative matrices, as described in the paper **ArXiv:2409.16822**. The algorithms provided here can be used to recreate all simulations presented in both the main paper and the supplementary material (or appendix, in the ArXiv version).

## Contents

- `illustrative_example.m`: A customizable script that applies all three algorithms to the same matrix family discussed in Section 5.1 of the paper, allowing for different values of M, target accuracy, maximum number of evaluations, etc.
- `real_antinorm.m`: Computes the polytope antinorm of a vector.
- `real_antinorm_cvx.m`: Computes the polytope antinorm of a vector relying on the CVX solver. See the instructions in the file since it requires installing CVX before using it.
- `aNorm.m`: Computes the polytope antinorm of a matrix. It can use real_antinorm (by default) or real_antinorm_cvx depending on the needs of the user.
- `lsr_comp_standard.m`: The standard algorithm, in which the polytope antinorm is fixed at the beginning and does not change.
- `adaptive_subradius_comp.m`: An adaptive variant that adds new vertices only if they improve the current antinorm.
- `adaptive_eigenvectors_subradius_comp.m`: Another adaptive variant that further refines the antinorm by adding suitably rescaled eigenvectors. Note: As mentioned in the paper, this is the only variant without a theoretical result guaranteeing convergence.

## Usage

The script `illustrative_example.m` can be run as it is and, if necessary, the customizable options can be easily modified following the comments in that files.

To apply the algorithms to a different matrix family, modify the `illustrative_example.m` file:

1. Replace `A_1` and `A_2` with the desired matrices.
2. Since the LSR is likely not known, renormalize the family using, for example, the result of `min(a(A_1),a(A_2))`, where `a()` is the chosen initial polytope antinorm. Alternatively, rescale the family using a preliminary bound which can be obtained by running `lsr_comp.m` with a maximum number of antinorm evaluations `M` small (e.g., equal to 10).
3. Modify `V_in` to define the new antinorm. For a family of dxd matrices, `V_in` should be a dxp matrix, where the p columns are the vertices. For example, `V_in = eye(d)` is the 1-antinorm in the dxd setting.

To apply the perturbation theory to a matrix family of dxd matrices:

1. Introduce the perturbation matrix `Delta = rand(d,d * m)`, were m is the number of matrices in the family (e.g., m = 2 above). Next, rescale the matrix by setting `Delta = Delta/norm(Delta,'fro')`.
2. Let `eps` be the size of the perturbation (e.g., between 10e-3 and 10e-7) and define the perturbed family `A_p = A + eps * Delta`.
3. Apply the adaptive algorithms following `illustrative_example.m` depending on the context to the perturbed family A_p. Possibly, do so for several values of `eps` for a more complete picture.
4. **Optional step:** Verify if the s.l.p. found for the perturbed family when `eps` is small is also a s.l.p. for the initial family when `eps = 0` using, for example, the polytopic algorithm.

## License

This project is licensed under the GNU General Public License v3.0. See the LICENSE file for more details.
