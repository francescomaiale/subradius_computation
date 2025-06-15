function [a, v] = aNorm(Y, V)

% aNorm.m Compute the polytope antinorm (defined by the vertex set V) of a matrix Y
%
%   [a, v] = ANORM(Y, V) computes, for each column V(:,i), the polytope antinorm
%   of the vector z = Y*V(:,i), then returns the smallest such bound across all i
%   and the corresponding transformed vector.
%
%   NOTE: This function calls real_antinorm.m, which must be on the MATLAB path.
%   If you prefer to use a CVX-based solver, replace real_antinorm with
%   real_antinorm_cvx (requires CVX toolbox: https://cvxr.com/cvx/), as 
%   explained in the corresponding MATLAB file real_antinorm_cvx.m
%
%   Inputs:
%     Y - The matrix to compute the antinorm of. Must satisfy size(Y,2)==size(V,1).
%     V - The matrix which columns define the polytope antinorm.
%
%   Outputs:
%     a - scalar equal to min_i lower_i, where lower_i is the lower bound
%         on the antinorm computed by real_antinorm(V, Y*V(:,i)).
%     v - the transformed vector z = Y*V(:,I) corresponding to the minimal bound.

    % Validate inputs
    if nargin < 2
        error('aNorm requires two inputs: Y (p x n) and V (n x j).');
    end
    [n, j] = size(V);
    [p, nY] = size(Y);
    if nY ~= n
        error('Inner dimensions must agree: size(Y,2) == size(V,1).');
    end

    tol = 1e-10;           % numerical tolerance for lower bound
    c = zeros(j, 1);       % preallocate for antinorm lower bounds

    % Loop over each column of V
    for i = 1:j
        z = Y * V(:, i);    % generate vector for real_antinorm
        try
            % Compute lower bound on real antinorm
            [lower_i, ~, ~] = real_antinorm(V, z);
            % If desired, swap to CVX solver:
            % [lower_i, ~, ~] = real_antinorm_cvx(V, z);
            c(i) = lower_i;
        catch ME
            warning('aNorm:SolverError', ...
                'Iteration %d: real_antinorm failed: %s. Setting bound = Inf.', i, ME.message);
            c(i) = Inf;
        end
    end

    % Find minimal bound and corresponding vector
    [a, idx] = min(c);
    v = Y * V(:, idx);
end
