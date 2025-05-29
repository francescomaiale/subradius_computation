%=============================================================================================================================================================
% Description:
%   Illustrative Example (Section 5.1) – compute the lower spectral radius
%   (LSR) of the matrix family F={A1,A2} via three methods:
%     • Standard Algorithm (S) (lsr_comp.m)
%     • Adaptive Algorithm (A) (adaptive_subradius_comp.m)
%     • Eigenvector-based Algorithm (E) (adaptive_eigenvectors_subradius_comp.m)
%
% Outputs:
%   • lsr_standard, p_standard (LSR bounds and performance metrics in Alg. S)
%   • lsr_adaptive_A, p_adaptive_A, n_vert_A (LSR bounds, performance metrics, and number of vertices at termination of the adaptive antinorm in Alg. A)
%   • lsr_adaptive_E, p_adaptive_E, n_vert_E (LSR bounds and performance metrics, and number of vertices at termination of the adaptive antinorm in Alg. E)
%=============================================================================================================================================================

%% 1. Define matrix family F = {A1, A2} and define A (2-by-4) to store each block [A1, A2] for convenience (passed to algorithms)
A_1 = [7 2; 0 3]; A_2 = [2 0; 4 8];
A = [A_1 A_2];

%% 2. Compute the LSR explicitly using the s.l.p. (Pi) through the formula lsr = spectral radius(Pi)^(1 / degree(Pi))
Pi = A_1*A_2*(A_1^2*A_2)^2;
rho = eigs(Pi,1,'largestabs')^(1/8);

%% 3. Renormalize family to avoid numerical overflow
A = A/rho;

%% 4. Set the common parameters for all algorithms
delta = 1e-6; % desired accuracy for LSR approximation
numeric_error = [0 0]; % assume entries of A are exact
display = 1;% show intermediate progress (set to false to suppress)
V_in = eye(2); % initial polytope antinorm set to the 1-antinorm

%% 5. Application of Algorithm S for different values of maximum allowed antinorm evaluations (M)
M_S = [500, 5e3, 1e4]; 
lsr_standard = zeros(numel(M_S),1);
p_standard   = cell(numel(M_S),1);

for idx = 1:numel(M_S)
    [lsr_standard(idx), p_standard{idx}, ~] = lsr_comp(A, numeric_error, delta, M_S(idx), V_in, [], display);
end

%% 6. Algorithm A (adaptive subradius computation with optional pruning)
%   - pruning flag toggles vertex pruning in the adaptive procedure
prune_A = true;                % set to false to disable pruning
M_A = [100 500 1e3];
lsr_adaptive_A = zeros(numel(M_A),1);
p_adaptive_A   = cell(numel(M_A),1);
n_vert_A       = zeros(numel(M_A),1);

for idx = 1:numel(M_A)
    [lsr_adaptive_A(idx), p_adaptive_A{idx}, ~, V_new] = adaptive_subradius_comp( A, numeric_error, delta, M_A(idx), V_in, [], display, prune_A);
    n_vert_A(idx) = size(V_new,2);  % number of vertices in the final polytope antinorm
end

%% 7. Algorithm E (adaptive eigenvector‐based subradius with optional pruning)
%   - theta: scaling parameter for eigenvectors prior to inclusion
%   - pruning flag toggles vertex pruning in the adaptive eigenvector method
theta   = 1.1;
prune_E = true;                % set to false to disable pruning
M_E         = [100, 500, 1e3];
lsr_adaptive_E = zeros(numel(M_E),1);
p_adaptive_E   = cell(numel(M_E),1);
n_vert_E       = zeros(numel(M_E),1);

for idx = 1:numel(M_E)
    [lsr_adaptive_E(idx), p_adaptive_E{idx}, ~, V_new] = adaptive_eigenvectors_subradius_comp( A, numeric_error, delta, M_E(idx), V_in, [], display, theta, prune_E);
    n_vert_E(idx) = size(V_new,2);
end

%% 8. Display summary of results
fprintf('\nResults summary (LSR estimates):\n');
fprintf('  Algorithm S:    M = [%g, %g, %g] → LSR ≈ [%g, %g, %g]\n', M_S, lsr_standard);
fprintf('  Algorithm A:    M = [%g, %g, %g] → LSR ≈ [%g, %g, %g];  #verts = [%d, %d, %d]\n', M_A, lsr_adaptive_A, n_vert_A);
fprintf('  Algorithm E:    M = [%g, %g, %g] → LSR ≈ [%g, %g, %g];  #verts = [%d, %d, %d]\n\n', M_E, lsr_adaptive_E, n_vert_E);
