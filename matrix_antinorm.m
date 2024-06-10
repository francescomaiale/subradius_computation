function [anorm,z] = matrix_antinorm(A,V)

% This function (Algorithm SM.1.2 in the supplementary material) computes the polytope antinorm (corresponding to the vertex set V) of the matrix A.
% The functon real_antinorm (corresponding to Algorithm SM.1.1) is called here to find the polytope antinorm of vectors.


%% Input
% A is a dxd real matrix for which we want to compute the antinorm
% V is a dxp real matrix which contains the p vertices defining the polytope antinorm as columns.


%% Output
% anorm is the antinorm of the matrix A
% z is the candidate new vertex in Algorithm A (lsr_comp_adaptive.m) and Algorithm E (lsr_comp_adaptive_eigenvectors.m) and coincides with A*v_i, where v_i is the column of V achieving the minimum in the computation of the antinorm


    
p = size(V,2); % Number of vertices
c = zeros(j,1); % Auxiliary vector to store the antinorm of A*v_i for all i = 1,...,p

for i = 1:p
    
    try
        
    [~,upper,~] = real_antinorm(V,Y*V(:,i)); % Antinorm of A*v_i
    c(i,1) = upper; % Store the value in the auxiliary vector
    
    catch ME % If there is any error, print a message and skip to the next iteration
        
        fprintf('Error on iteration %d: %s\n', i, ME.message);
        continue;
        
    end
    
end

[anorm,I] = min(c(:,1)); % Compute the minimum value of a(A*v_i) with respect to i and denote by I the index such that min = a(A*v_I)

z = Y*V(:,I); % Candidate new vertex

end