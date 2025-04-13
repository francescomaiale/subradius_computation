function [a, v] = aNorm(Y, V)
    % aNorm  Compute the "antinorm" of a matrix V under transformation Y.
    %
    % Inputs:
    %   Y : an (n x n) or (n x k) matrix
    %   V : an (n x j) matrix whose columns are vectors for which the antinorm is computed
    %
    % Outputs:
    %   a : the minimum "upper" value of real_antinorm(V, Y*V(:,i)) over all columns i
    %   v : the vector Y * V(:,I) that achieves this minimum

    j = size(V, 2);      % number of columns (vectors) in V
    c = zeros(j, 1);     % preallocate array for the "upper" values

    for i = 1 : j
        try
            [lower, ~, ~] = real_antinorm(V, Y * V(:, i));
            c(i) = lower;
        catch ME
            % If real_antinorm fails for any reason, you can choose to skip or handle differently
            fprintf('Error on iteration %d: %s\n', i, ME.message);
            c(i) = Inf; 
            continue;
        end
    end

    [a, I] = min(c);

    v = Y * V(:, I);

end
