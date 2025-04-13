function [lower, upper, x] = real_antinorm(V,z)

    [n,m] = size(V);

    A_ineq = [ -z,      V;
                0, -ones(1,m) ];
    b_ineq = [ zeros(n,1);
               -1 ];

    c = zeros(m+1,1);
    c(1) = 1;

    lb = zeros(m+1,1);
    ub = [];

    Aeq = [];
    beq = [];
    tol = 1e-10;  

    options = optimset('Algorithm','dual-simplex-highs','MaxIter',max(300,n), 'Display','off','TolFun',tol);
    [x, ~, exitflag] = linprog(c, A_ineq, b_ineq, Aeq, beq, lb, ub, options);

    lower = 0;
    upper = 0;

    switch exitflag

        case 1
            % ------------------------------------------------------------------
            % *Optimal feasible solution found*
            % ------------------------------------------------------------------
            if isempty(x)
                return
            end

            xVal = x(1);

            if xVal <= tol
                lower = Inf;
                upper = Inf;
                return
            end

            denomLower = xVal + tol;
            denomUpper = xVal - tol;

            if denomUpper <= 0
                lower = Inf;
                upper = Inf;
            else
                upper = 1 / denomUpper;
                lower = 1 / denomLower;
            end

        case 0
            % ------------------------------------------------------------------
            % *Reached iteration or time limit*
            % ------------------------------------------------------------------
            warning('linprog:NoConvergence',...
                    'Iteration limit or time limit reached, solution may be inaccurate.');
            if ~isempty(x) && x(1) > 0
                denomLower = x(1) + tol;
                denomUpper = x(1) - tol;
                if denomUpper <= 0
                    lower = Inf; 
                    upper = Inf;
                else
                    upper = 1 / denomUpper;
                    lower = 1 / denomLower;
                end
            else
                lower = 0; 
                upper = 0;
            end

        case -2
            % ------------------------------------------------------------------
            % *No feasible solution* (infeasible)
            % ------------------------------------------------------------------

            lower = 0;
            upper = 0;
            return

        case -3
            % ------------------------------------------------------------------
            % *Problem is unbounded*
            % ------------------------------------------------------------------
            lower = 0;  % or 0
            upper = Inf;
            return

        otherwise
            % ------------------------------------------------------------------
            % *Other negative or unexpected EXITFLAGs*
            % ------------------------------------------------------------------
            warning('linprog:Exitflag',...
                'linprog returned EXITFLAG = %d, solution may be invalid.', exitflag);
            lower = 0;
            upper = 0;
            return
    end
end
