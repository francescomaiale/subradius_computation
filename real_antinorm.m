function [lower, upper, x] = real_antinorm(V, z)
% real_antinorm.m Compute bounds on the polytopic antinorm (defined by the vertices
% given by the columns of V) of a vector z.
%
%   [lower, upper, x] = real_antinorm(V, z) solves the following linear program:
%       minimize    t
%       subject to  V * alpha <= t * z,
%                   sum(alpha) >= 1,
%                   alpha >= 0.
%
%   Inputs:
%     V      - n-by-m real matrix.
%     z      - n-by-1 real vector.
%
%   Outputs:
%     lower  - lower bound on the antinorm (scalar). If the problem is
%              unbounded or degenerates, this may be Inf or 0 as indicated.
%     upper  - upper bound on the antinorm (scalar). If the problem is
%              unbounded or degenerates, this may be Inf.
%     x      - LP solution vector of length m+1, structured as [t; alpha],
%              where t = x(1) and alpha = x(2:end).

    % Validate inputs
    if nargin < 2
        error('real_antinorm requires two inputs: V (vertex matrix) and z (vector).');
    end
    [n, m] = size(V);
    if ~isvector(z) || length(z) ~= n
        error('z must be a column vector of length matching the rows of V.');
    end

    % Set up inequality constraints for linprog:
    %   -z*t + V*alpha <= 0  corresponds to  [-z, V] * [t; alpha] <= 0
    %   sum(alpha) >= 1      corresponds to  [0, -ones(1,m)] * [t; alpha] <= -1
    A_ineq = [ -z,      V;
                 0, -ones(1,m) ];
    b_ineq = [ zeros(n,1);
               -1 ];

    % Minimize t (scalar), so c = [1; zeros(m,1)]
    c = zeros(m+1,1);
    c(1) = 1;

    % Variable bounds: t >= 0, alpha >= 0
    lb = zeros(m+1,1);
    ub = [];

    % No equality constraints
    Aeq = [];
    beq = [];

    % Numerical tolerance
    tol = 1e-10;

    % Set solver options (prefer Highs dual-simplex, else fallback)
    try
        options = optimset('Algorithm','dual-simplex-highs', ...
                           'MaxIter', max(300, n), ...
                           'Display','off', ...
                           'TolFun', tol);
    catch
        options = optimset('Algorithm','dual-simplex', ...
                           'MaxIter', max(300, n), ...
                           'Display','off', ...
                           'TolFun', tol);
    end

    % Solve the LP
    [x, ~, exitflag] = linprog(c, A_ineq, b_ineq, Aeq, beq, lb, ub, options);

    % Initialize outputs
    lower = 0;
    upper = 0;

    % Interpret solver exit flags
    switch exitflag
        case 1  % Optimal solution found
            if isempty(x)
                return
            end
            t_val = x(1);
            if t_val <= tol
                lower = Inf;
                upper = Inf;
                return
            end
            % Compute reciprocal bounds with safety tol adjustment
            denomLow  = t_val + tol;
            denomHigh = t_val - tol;
            if denomHigh <= 0
                lower = Inf;
                upper = Inf;
            else
                upper = 1 / denomHigh;
                lower = 1 / denomLow;
            end

        case 0  % Iteration or time limit reached
            warning('linprog:NoConvergence', ...
                    'Iteration or time limit reached; solution may be inaccurate.');
            if ~isempty(x) && x(1) > tol
                denomLow  = x(1) + tol;
                denomHigh = x(1) - tol;
                if denomHigh <= 0
                    lower = Inf;
                    upper = Inf;
                else
                    upper = 1 / denomHigh;
                    lower = 1 / denomLow;
                end
            else
                lower = 0;
                upper = 0;
            end

        case -2  % No feasible solution
            lower = 0;
            upper = 0;
            return;

        case -3  % Problem is unbounded
            lower = 0;
            upper = Inf;
            return;

        otherwise  % Other unexpected flags
            warning('linprog:Exitflag', ...
                    'linprog returned exitflag = %d; solution may be invalid.', exitflag);
            lower = 0;
            upper = 0;
            return;
    end
end
