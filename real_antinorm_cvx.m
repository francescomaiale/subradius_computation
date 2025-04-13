function [lower, upper, x_sol] = real_antinorm_cvx(V, z)

    [n, m] = size(V);

    tol = 1e-10;

    lower  = 0;
    upper  = 0;
    x_sol  = [];

    %----------------------- CVX PROBLEM SETUP -----------------------%
    cvx_begin
        cvx_precision best
        
        variable x_cvx(m+1)  % x_cvx(1) = the scalar part, x_cvx(2:end) = combination coeffs

        minimize( x_cvx(1) )

        for i = 1:n
            -z(i)*x_cvx(1) + V(i,:)*x_cvx(2:end) <= 0;
        end
        -sum(x_cvx(2:end)) <= -1;  

        x_cvx >= 0;
    cvx_end

    x_sol = x_cvx;

    status = cvx_status;
    
    switch status
        case 'Solved'
            % --------------------------------------------------------------
            % *Optimal, feasible solution found*
            % --------------------------------------------------------------
            if isempty(x_cvx)
                return
            end
            xval = x_cvx(1);

            if xval <= tol
                lower = Inf;
                upper = Inf;
                return
            end

            denom_upper = xval - tol;  
            denom_lower = xval + tol; 

            if denom_upper <= 0
                lower = Inf;
                upper = Inf;
            else
                upper = 1 / denom_upper;
                lower = 1 / denom_lower;
            end

        case 'Infeasible'
            % --------------------------------------------------------------
            % *No feasible solution*
            % --------------------------------------------------------------
            lower = 0;
            upper = 0;

        case 'Unbounded'
            % --------------------------------------------------------------
            % *Problem is unbounded*
            % --------------------------------------------------------------
            lower = 0;   % or 0
            upper = Inf;

        case {'Inaccurate', 'Inaccurate/Solved'}
            % --------------------------------------------------------------
            % *Solution found but CVX indicates potential numerical issues*
            % --------------------------------------------------------------
            warning('real_antinorm_cvx:Inaccurate',...
                'CVX reports solution may be inaccurate. Proceeding with caution.');
            if ~isempty(x_cvx) && x_cvx(1) > tol
                denom_upper = x_cvx(1) - tol;
                denom_lower = x_cvx(1) + tol;
                if denom_upper <= 0
                    lower = Inf;
                    upper = Inf;
                else
                    upper = 1 / denom_upper;
                    lower = 1 / denom_lower;
                end
            else
                lower = 0;
                upper = 0;
            end

        otherwise
            % --------------------------------------------------------------
            % *Failed or other unexpected solver status*
            % --------------------------------------------------------------
            warning('real_antinorm_cvx:Failed',...
                'CVX returned status ''%s''; solution may be invalid.', cvx_status);
            lower = 0;
            upper = 0;
    end
end
