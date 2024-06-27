function [lsr,p,err,vertex_set] = adaptive_subradius_comp(A,dA,delta,M,vertex_set,err,display)
% This file is part of the subradius_computation project.
% Copyright (c) 2024
% This source code is licensed under the GNU General Public License v3.0. See the LICENSE file for more details.


% This function (Algorithm (A) in the main paper) computes a lower and upper bound for the 
% lower spectral radius of a matrix family {A_1, ..., A_m} starting from an initial polytope antinorm and 
% refining it adaptively.


%% Input
% A = [A_1, ..., A_m], given by horizontal concatenation of the matrices
% dA is an estimate of the norm of the error in the matrices in A
% delta = is the pre-selected precision within which the LSR should be estimated; specifically, (upper bound - lower bound) should be less than delta
% M is the maximum number of antinorm evalutions
% vertex_set is the vertex set corresponding to the initial polytope antinorm
% err = [err(1), err(2)] gives the relative errors in evaluating products and antinorms
% display is either 1 or 0; if 1, the algorithm prints the performance vector p and the lower and upper bound after each iteration of the main while loop

%% Output
% lsr = [lower bound for the LSR, upper bound for the LSR]
% p is the performance metrics vector, namely p = [degree for which the gap between upper and lower bound is optimal, degree which yields the best upper bound, number of evaluations of antinorms, maximum number of matrices considerd at each stage]
% vertex_set is the final set of vertices, which corresponds to the refined polytope antinorm at when the algorithm ends


%% Check if there are any errors in the input values   
if nargin < 2
    dA = [];
end

if nargin < 3
    delta = 1e-6;
end

if nargin <4
    M = 500;
end

% If no initial antinorm is selected, we use the 1-antinorm
if nargin < 5
    vertex_set = eye(size(A,1));
end

if nargin < 6
    err = [];
end

if nargin < 7
    display = 0;
end


%% Check if the matrix A has the correct dimension
[ma,na] = size(A); % ma = number of rows (or columns) of each A_i, so if the input is correct na/ma should coincide with the number of elements of the family

if (na > ma && rem(na,ma) == 0)
    m = fix(na/ma); % number of elements of the family
    d = ma; % each matrix in the family is (d x d)
elseif (ma > na && rem(ma,na) == 0) % happens when the input A has been given as vertical concatenation of the matrices A_i. Considering A transpose solves the issue:
    A = A';
    m = fix(ma/na);
    d = na;
else % if neither situation happens, there is an error in the input A
    disp('The matrix has wrong dimensions');
    return;
end

%% Adjust the relative errors in evaluating products and antinorms and the error in the norm of A
if  size(err)*[1;1] < 3
    err = [2^(-50)*d,2^(-49)*d];
end

dA = dA(:);

if  size(dA)*[1;0] == d
    proportional = 0;
else 
    proportional = 1;
    if size(dA)*[1;0] == 0
        dA = eps;
    end
end


%% First step: evaluation of products of length one (i.e., the matrices A_i)
tol = 1e-10; % tolerance used for 1) criterion to include new vertices, 2) pruning of the vertex set.

X = []; % initialize to store matrices
H_Radius = inf; % initial guess for the upper bound
L_Radius = 0; % initial guess for the lower bound

antinorm = zeros(m,1); lowbound = zeros(m,1); antinorm_err = zeros(m,1); % preallocate space

for i = 1:m
    Y = A(:,(i-1)*d+1:i*d); % set Y := A_i
    [aN,v_new]=matrix_antinorm(Y,vertex_set); % compute the antinorm of Y with respect to the polytope antinorm which correponds to the vertex set; 
                                              % v_new, on the other hand, is the candidate vertex that solves a(Y)= min_{v in vertex_set} a(Yv)
  
    antinorm(i) = aN*(1+err(2)); 
    lowbound(i) = antinorm(i);
  
    if proportional
        antinorm_err(i) = dA(1) * antinorm(i);
    else
        antinorm_err(i) = dA(i);
    end
    
    H_Radius = min(H_Radius,abs(eigs(Y,1,'largestabs'))); % update the upper bound as the minimum between the previous upper bound and the spectral radius of Y=A_i

    % Adaptive refinement of the vertex set/antinorm
    [~,upper,~] = real_antinorm(vertex_set,v_new); % Compute the antinorm of the candidate new vertex, v_new, with respect to the current vertex set
    % If the antinorm of v_new is less than or equal to one, it falls inside (or on the boundary) of the polytope; thus, we add it to the vertex set. We use a small tolerance to take into account possible numerical inaccuracies
    if upper <= 1+tol
        vertex_set = [vertex_set v_new];
    end
    
    X = [X;Y]; % add Y to the storage X
end
 

%% Vertices pruning (Algorithm SM.1.3 in the supplementary material)

reduce=0; % parameter used to exit the while loop if no vertex is removed in an entire for loop

while (reduce == 0)
    n_v = size(vertex_set,2); % current number of vertices
    
    if (n_v <= 1) % exit if there is only one vertex left
        break;
    end
    
    reduce = 1;
    
    for i = n_v:-1:1
        v = vertex_set(:,i); % consider the i-th vertex, i.e. i-th column of the vertex set
        W = [vertex_set(:,1:i-1) vertex_set(:,i+1:n_v)]; % construct the auxiliary vertex set, obtained from the 
        if (rank(W)==1)
            break
        end
        [lower,~,~] = real_antinorm(W,v);
        if (lower >= 1 + tol)
            reduce=0;
            vertex_set = W;           
            break;        
        end        
    end
end
% Find indices of unique columns and extract them
[~, idx] = unique(vertex_set', 'rows', 'stable'); vertex_set = vertex_set(:, idx);


%% Initialize the iteration parameter and update the lower bound:
dA = antinorm_err;
antinorm_A = antinorm;
L_Radius = max(L_Radius,min(antinorm+antinorm_err));

n_op = m; % number of antinorm evaluations a(P) made so far corresponds with the number of matrices in the family.
n = 1; % current degree of matrices that have been evaluated by the algorithm
J = m; % number of elements evaluated among products of degree n. Since n = 1, it corresponds to the number of matrices in the family.
MaxJ = J; % keeps track of the maximum value achieved by J

BestPower = 1; % keeps track of the degree yielding the optimal gap between lower and upper bound
ell_slp = 1; % keeps track of the degree yielding the optimal upper bound

%% Main loop of the algorithm
while n_op < M && L_Radius < H_Radius - delta
    n = n + 1; % increase n as we start to explore products of degree n + 1
    H_Radius_Old = H_Radius; L_Radius_Old = L_Radius; % store the previous upper and lower bound to later establish if 1) the gap has improved and 2) the upper bound has improved
    L_Radius = H_Radius; 
    NewJ = 0; % this variable counts how many products of degree n + 1 are evaluated by the algorithm

    % The first "for" cycle considers products of length n-1, evaluated already in the previous step. 
    % The second "for" cycle considers the elements of the family, which multiplied with products of length n-1 generate products of length n

    for k = 1:J
        for i = 1:m
            NJJN = NewJ + 1;
            Y = X((k-1)*d+1:k*d,:)*A(:,(i-1)*d+1:i*d); % current product of length n
            
            antinorm_error_new(NJJN) = (antinorm(k) * (err(1) * antinorm_A(i) + dA(i)) + antinorm_err(k) * (antinorm_A(i) + dA(i)));
            [aN,v_new] = matrix_antinorm(Y, vertex_set); % compute the antinorm of Y
            antinorm_new(NJJN) = aN*(1+err(2));
            lowbound_new(NJJN) = max(lowbound(k),(antinorm_new(NJJN)+antinorm_error_new(NJJN))^(1/n)); % update the lowbound vector
            % Add the candidate vertex if it satisfies the criterion (as in the first step)
            [~,upper,~] = real_antinorm(vertex_set,v_new);
            if upper <= 1+tol
                vertex_set = [vertex_set v_new];
            end

            H_Radius = min(H_Radius, (abs(eigs(Y,1,'largestabs')))^(1/n)); % update the upper bound as the min between the current value and the spectral radius
            if lowbound_new(NJJN) < H_Radius - delta % establish if the matrix Y make the cut and will be used in the next step to generate products of degree n+1
                NewJ = NewJ + 1; 
                L_Radius = min(L_Radius,lowbound_new(NJJN)); % update the lower bound
                NewX((NewJ-1)*d+1:NewJ*d,:) = Y; % store this matrix as it will be used in the next iteration of the while loop
            end
        end
    end

    % If no matrix product satisfies the "if" above, then convergence to delta has been achieved; thus, terminate the algorithm
    if isempty(NewX)
        return
    else
    X = NewX; % replace the matrices stored previosly (in X) with the new ones to prepare for the next iteration.
    end

    L_Radius = max(L_Radius_Old,min(L_Radius,H_Radius-delta)); % final update of the lower bound
    n_op = n_op + J * m; % each product yields one antinorm evaluation, so the current number increases by the product of the lengths of the two for cycles.
    
    % Replace old antinorms, lowbounds and errors with the new ones
    antinorm = antinorm_new; lowbound = lowbound_new; antinorm_err = antinorm_error_new;
    
    % Prepare J for the next iteration and update MaxJ
    J = NewJ; MaxJ = max(J,MaxJ);

    if H_Radius - L_Radius < H_Radius_Old - L_Radius_Old
        BestPower = n; % update the optimal gap degree if the gap decreased
    end

    if H_Radius < H_Radius_Old
        ell_slp=n; % update the optimal upper bound degree if the upper bound decreased
    end
  
    if display
        disp([BestPower,ell_slp,n,n_op,J]);
        disp([L_Radius,H_Radius]);
    end

%% Pruning of vertices before the next iteration of the main loop
    reduce=0; % Initialize the pruning
    while (reduce==0)
        n_v =size(vertex_set,2);
        if (size(vertex_set,2)<=1)
            break;
        end
        reduce=1;
        for l=n_v:-1:1
            v=vertex_set(:,l);
            W=[vertex_set(:,1:l-1) vertex_set(:,l+1:n_v)];
            if (rank(W)==1)
                break
            end
            [lower,~,~] = real_antinorm(W,v); % In this instance, the lower bound only is needed
            if (lower>=1+tol)
                %disp('eliminate vertex');
                reduce=0;
                vertex_set=W;             
                break;        
            end        
        end
    end
[~, idx] = unique(vertex_set', 'rows', 'stable'); vertex_set = vertex_set(:, idx);
end

p = [BestPower,ell_slp,n,n_op,MaxJ];
lsr = [L_Radius,H_Radius];
