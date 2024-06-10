function [anorm,z] = matrix_antinorm(A,V)

% This function (Algorithm SM.1.2 in the supplementary material) computes the polytope antinorm (corresponding to the vertex set V) of the matrix A.
% The functon real_antinorm (corresponding to Algorithm SM.1.1) is called here to find the polytope antinorm of vectors.


%% Input
% A is a dxd real matrix for which we want to compute the antinorm
% V is a dxp real matrix which contains the p vertices defining the polytope antinorm as columns


%% Output
% anorm = antinorm of the matrix A
% z is the candidate new vertex in Algorithm A (adaptive_subradius_comp.m) and Algorithm E (adaptive_eigenvectors_subradius_comp.m)
% and coincides with A*v_i, where v_i is the column of V achieving the minimum in the computation of the antinorm


    
p = size(V,2); % number of vertices
c = zeros(p,1); % auxiliary vector that stores the antinorm of A*v_i for all i = 1,...,p

for i = 1:p
    
    try
        
    [~,upper,~] = real_antinorm(V,A*V(:,i)); % compute the antinorm of A*v_i
    c(i,1) = upper; % store the value in the auxiliary vector
    
    catch ME % if there is any error, print a message and skip to the next iteration
        
        fprintf('Error on iteration %d: %s\n', i, ME.message);
        continue;
        
    end
    
end

[anorm,I] = min(c(:,1)); % compute the minimum value of a(A*v_i) with respect to i and denote by I the index such that min = a(A*v_I)

z = A*V(:,I); % candidate new vertex

end
