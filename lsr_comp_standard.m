function [lsr,p,err] = lsr_comp_standard(A, dA, delta, M, V, err, display)

% This function computes a lower and upper bound for the lower spectral radius of a matrix family {A_1, ..., A_m}.
% It corresponds to Algoritm 4.1 in the main paper.

%% Input
% A = [A_1, ..., A_m], given by horizontal concatenation of the matrices
% dA is an estimate of the norm of the error in the matrices in A
% delta = is the pre-selected precision within which the LSR should be estimated; specifically, (upper bound - lower bound) should be less than delta
% M is the maximum number of antinorm evalutions
% V is the vertex set of the polytope antinorm
% err = [err(1), err(2)] gives the relative errors in evaluating products and antinorms

%% Output
% lsr = [lower bound for the LSR, upper bound for the LSR]
% p is the performance metrics vector, namely p = [degree for which the gap between upper and lower bound is optimal, degree which yields the best upper bound, number of evaluations of antinorms, maximum number of matrices considerd at each stage]
    

%% Check if there are any mistakes in the input values

if nargin < 2
  dA = [];
end

if nargin < 3
  delta = 1e-6 
end

if nargin < 4
  M = 500;
end

% If no vertex set is given to this function, it utilizes the 1-antinorm (which vertex set is given by (1,0,...,0), ..., (0,...,0,1))
if nargin < 5
    V = eye(size(A,1));
end

if nargin < 6
  err = [];
end

if nargin < 7
    display=1;
end

%% Check if the matrix A has the correct dimension

[ma,na] = size(A); % ma = number of rows (or columns) of each A_i, so if the input is correct na/ma should coincide with the number of elements of the family

if (na > ma && rem(na,ma) == 0)
  m = fix(na/ma); % Number of elements of the family
  d = ma; % Each matrix in the family is (d x d)
  
elseif (ma > na && rem(ma,na) == 0) % This happens if the input A has been given as vertical concatenation of the matrices A_i. Considering A transpose solves the issue:
  A = A';
  m = fix(ma/na);
  d = na;
  
else % If neither situation happens, there is an error in the input A.
 disp('The matrix has wrong dimensions');
 return;
end

%% Adjust the relative errors in evaluating products and antinorms and the error in the norm of A.

if  size(err)*[1;1] < 3
    err = [2^(-50)*d,2^(-49)*d];
end

dA = dA(:);

if  size(dA)*[1;0] == d
  prop = 0;
else 
  prop = 1;
  if size(dA)*[1;0] == 0 % If the input is dA = [], then dA is set equal to the eps-precision of the machine.
    dA = eps;
  end
end


%% First step: evaluation of products of length one (i.e., the matrices A_i)

X = []; % Initialize to store matrices

H_Radius = inf; % Initial guess for the upper bound
L_Radius = 0; % Initial guess for the lower bound

antinorm = zeros(m,1); lowbound = zeros(m,1); antinorm_err = zeros(m,1); % Preallocate space

for i = 1:m

  Y = A(:,(i-1)*d+1:i*d); % Set Y = A_i
  
  [aN,~] = matrix_antinorm(Y,V); % Compute the antinorm of Y with respect to the polytope antinorm which correponds to the vertex set V
  
  antinorm(i) = aN*(1+err(2)); 
  lowbound(i) = antinorm(i);
  
  if prop
    antinorm_err(i) = dA(1) * antinorm(i);
  else
    antinorm_err(i) = dA(i);
  end

  H_Radius = min(H_Radius,abs(eigs(Y,1,'largestabs'))); % Update the upper bound as the minimum between the previous upper bound and the spectral radius of Y=A_i
 
  X = [X;Y]; % Add Y to the storage X

end


%% Initialize the iteration parameter and update the lower bound:

dA = antinorm_err;
antinorm_A = antinorm;
L_Radius = max(L_Radius,min(antinorm+antinorm_err));

n_op = m; % Number of antinorm evaluations a(P) made so far corresponds with the number of matrices in the family.
n = 1; % Current degree of matrices that have been evaluated by the algorithm
J = m; % Number of elements evaluated among products of degree n. Since n = 1, it corresponds to the number of matrices in the family.
MaxJ = J; % Keeps track of the maximum value achieved by J

BestPower = 1; % Keeps track of the degree yielding the optimal gap between lower and upper bound
ell_slp = 1; % Keeps track of the degree yielding the optimal upper bound


%% Main loop of the algorithm

while n_op < M && L_Radius < H_Radius - delta

  n = n + 1; % Increase n as we start to explore products of degree n + 1

  H_Radius_Old = H_Radius; L_Radius_Old = L_Radius; % Store the previous upper and lower bound to later establish if 1) the gap has improved and 2) the upper bound has improved.
  
  L_Radius = H_Radius; 

  NewJ = 0; % Use this variable to count how many products of degree n + 1 are evaluated by the algorithm.

 % The first "for" cycle considers products of length n-1, evaluated already in the previous step. The second "for" cycle considers the elements of the family, which multiplied with products of length n-1 generate products of length n.
 
  for k = 1:J
    for i = 1:m

     NJJN = NewJ + 1; % This keeps track of all generated product of length n. It is different from NewJ because NewJ counts product that make the cut of the algorithm and will be used in the next step to generate products of length n+1; NJJN, on the other hand, counts all products of length n generated here, even those that are discarded later.
     
     Y = X((k-1)*d+1:k*d,:)*A(:,(i-1)*d+1:i*d); % Product of length n

     antinorm_error_new(NJJN) = (antinorm(k)*(err(1)*antinorm_A(i)+dA(i)) + antinorm_err(k)*(antinorm_A(i)+dA(i)));

     [aN,~] = matrix_antinorm(Y,V); % Compute the antinorm of Y
     antinorm_new(NJJN) = aN*(1+err(2));

     lowbound_new(NJJN) = max(lowbound(k),(antinorm_new(NJJN)+antinorm_error_new(NJJN))^(1/n)); % update the lowbound vector

     H_Radius = min(H_Radius, (abs(eigs(Y,1,'largestabs')))^(1/n)); % Update the upper bound as the min between the current value and the spectral radius

     if lowbound_new(NJJN) < H_Radius - delta % This decides if the matrix Y make the cut and will be used in the next step to generate products of degree n+1 or not.
         
        NewJ = NewJ + 1; 
        L_Radius = min(L_Radius,lowbound_new(NJJN)); % Update the lower bound
        NewX((NewJ-1)*d+1:NewJ*d,:) = Y; % Store this matrix as it will be used in the next iteration of the while loop
     end

    end
 end

  X = NewX; % Replace the matrices stored previosly (X) with the new ones to prepare for the next iteration.

  L_Radius = max(L_Radius_Old,min(L_Radius,H_Radius-delta)); % Update of the lower bound

  n_op = n_op + J*m; % Each product yields one antinorm evaluation, so the current number increases by the product of the lengths of the two for cycles.
  
  % replace old antinorms, lowbounds and errors with the new ones
  antinorm = antinorm_new; lowbound = lowbound_new; antinorm_err = antinorm_error_new;
    
  % prepare J for the next iteration and update MaxJ
  J = NewJ;
  MaxJ = max(J,MaxJ);

  if H_Radius - L_Radius < H_Radius_Old - L_Radius_Old
     BestPower = n; % Update the optimal gap degree if the gap decreases
  end

  if H_Radius < H_Radius_Old
      ell_slp=n; % Update the optimal upper bound degree if the upper bound decreases
  end
  
  if display
      disp([BestPower,ell_slp,n,n_op,J]);
      disp([L_Radius,H_Radius]);
  end
  
end

p = [BestPower,ell_slp,n,n_op,MaxJ];
lsr = [L_Radius,H_Radius];