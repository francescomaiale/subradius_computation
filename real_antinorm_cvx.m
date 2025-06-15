function [lower, upper, x_sol] = real_antinorm_cvx(V, z)

% real_antinorm_cvx.m Compute bounds on the polytope antinorm using CVX
% specifically, it solves the same antinorm LP as REAL_ANTINORM but relies on CVX.
%
%   Requires CVX toolbox: https://cvxr.com/cvx/
%
%   Inputs:
%     V      - n-by-m real matrix of vertices (columns).
%     z      - n-by-1 real vector.
%
%   Outputs:
%     lower  - lower bound on the antinorm (scalar). Can be Inf or 0.
%     upper  - upper bound on the antinorm (scalar). Can be Inf or 0.
%     x_sol  - solution vector of length m+1: [t; alpha], where t = x_sol(1)
%              and alpha = x_sol(2:end).
%
%   This formulation uses CVX to declare variables & constraints:
%       minimize   t
%       subject to V*alpha <= t*z,
%                  sum(alpha) >= 1,
%                  alpha >= 0.
%   Results are post-processed similarly to real_antinorm.m to ensure numeric safety.

    % Validate inputs
    if nargin < 2
        error('real_antinorm_cvx requires two inputs: V (vertices matrix) and z (vector).');
    end
    [n, m] = size(V);
    if ~isvector(z) || length(z) ~= n
        error('z must be a column vector matching the rows of V.');
    end

    % Numerical tolerance
    tol = 1e-10;

    % Initialize outputs
    lower = 0;
    upper = 0;
    x_sol = [];

    %----------------------- CVX PROBLEM SETUP -----------------------%
    % Ensure CVX is installed and on the path
    % https://cvxr.com/cvx/
    cvx_begin
        cvx_precision best          % maximize numeric accuracy, but it can be changed depending on the needs
        variable x_cvx(m+1)        % x_cvx(1)=t, x_cvx(2:end)=alpha

        % Objective: minimize scalar t
        minimize( x_cvx(1) )

        % Constraints:
        %   -z(i)*t + V(i,:)*alpha <= 0  <=> V(:,j)*alpha <= t*z
        for i = 1:n
            -z(i)*x_cvx(1) + V(i,:)*x_cvx(2:end) <= 0;
        end
        % Convex combination of columns must sum to at least 1
        -sum(x_cvx(2:end)) <= -1;

        % Nonnegativity of t and alpha
        x_cvx >= 0;
    cvx_end

    % Retrieve solution and status
    x_sol = x_cvx;
    status = cvx_status;

    % Interpret CVX status
    switch status
        case 'Solved'
            % Optimal solution found
            if isempty(x_cvx)
                return;
            end
            t_val = x_cvx(1);
            % infinite bounds
            if t_val <= tol
                lower = Inf;
                upper = Inf;
                return;
            end
            % Compute reciprocal bounds safely
            denom_lo = t_val + tol;
            denom_hi = t_val - tol;
            if denom_hi <= 0
                lower = Inf;
                upper = Inf;
            else
                upper = 1 / denom_hi;
                lower = 1 / denom_lo;
            end

        case 'Infeasible'
            % No feasible solution
            lower = 0;
            upper = 0;

        case 'Unbounded'
            % Problem is unbounded
            lower = 0;
            upper = Inf;

        case {'Inaccurate', 'Inaccurate/Solved'}
            % Numerical issues flagged
            warning('real_antinorm_cvx:Inaccurate', ...
                'CVX reports solution may be inaccurate.');
            if ~isempty(x_cvx) && x_cvx(1) > tol
                denom_hi = x_cvx(1) - tol;
                denom_lo = x_cvx(1) + tol;
                if denom_hi <= 0
                    lower = Inf;
                    upper = Inf;
                else
                    upper = 1 / denom_hi;
                    lower = 1 / denom_lo;
                end
            else
                lower = 0;
                upper = 0;
            end

        otherwise
            % Other failure or unexpected status
            warning('real_antinorm_cvx:Failed', ...
                'CVX returned status ''%s''; solution may be invalid.', status);
            lower = 0;
            upper = 0;
    end
end
