% In this file, we consider the example of Section 5.1 (Illustrative Example) and apply all algorithms: 
% Algorithm S (lsr_comp.m), Algorithm A (adaptive_subradius_comp.m), and Algorithm E (adaptive_eigenvectors_subradius_comp.m)

A_1 = [7 2; 0 3]; A_2 = [2 0; 4 8];
A = [A_1 A_2];

Pi = A_1*A_2*(A_1^2*A_2)^2; % known s.l.p. for the family 
rho = eigs(Pi,1,'largestabs')^(1/8); % lsr = spectral radius(s.l.p.)^(1 / degree(s.l.p.))

% For simplicity and to avoid large entries appearing when we let the algorithms investigate
% very long products, we apply our algorithms to the renormalized family:
A = A/rho; 

delta = 1e-6; % pre-selected accuracy to which we want to approximate the lsr
numeric_error = [0 0]; % since there is no numerical error on the entries of A
display = 1; % to show the progress of the algorithms after each step; set display=0 otherwise
V_in = eye(2); % 1-antinorm is the initial polytope antinorm in all cases


%% Algorithm S (Algorithm 4.1 in the main paper)
M=[500 5000 10000]; % different values of maximum number of allowed antinorm evaluations

for i = 1:length(M)
    [lsr_standard(i,:),p_standard(i,:),~] = lsr_comp(A, numeric_error, delta, M(i), V_in, [], display);
end

%% Algorithm A
M = [100 500 1000];
n_vert_A = zeros(length(M),1); % store the number of vertices of the adaptive antinorms
V_adaptive_A = []; % store the vertices of the adaptive antinorms

for i = 1:length(M)
    [lsr_adaptive_A(i,:),p_adaptive_A(i,:),~,V_new] = adaptive_subradius_comp(A, numeric_error, delta, M(i), V_in, [], display, true);
    V_adaptive_A = [V_new zeros(2,1)]; % the column (0, 0)^T is added to separate the different vertices set obtained for different values of M.
    n_vert_A(i) = size(V_new,2);
end


%% Algorithm E
M = [100 500 1000];
n_vert_E = zeros(length(M),1); % store the number of vertices of the adaptive antinorms
V_adaptive_E = []; % store the adaptive antinorms
theta = 1.1; % scaling parameter (to rescale eigenvectors prior to adding them)

for i = 1:length(M)
    [lsr_adaptive_E(i,:),p_adaptive_E(i,:),~,V_new] = adaptive_eigenvectors_subradius_comp(A, numeric_error, delta, M(i), V_in, [], display,theta,true);
    V_adaptive_E = [V_new zeros(2,1)];
    n_vert_E(i) = size(V_new,2);
end
