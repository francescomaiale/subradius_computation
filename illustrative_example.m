% In this file, we consider the example of Section 5.1 (Illustrative Example) and apply all our algorithms: Algorithm S (lsr_comp_standard.m), Algorithm A (adaptive_subradius_comp.m), and Algorithm E (adaptive_eigenvectors_subradius_comp.m)

close all
clearvars

%% Setting of the example

A_1 = [7 2; 0 3]; A_2 = [2 0; 4 8];

A = [A_1 A_2];

Pi = A_1*A_2*(A_1^2*A_2)^2; % s.l.p. for the family

rho = eigs(Pi,1,'largestabs')^(1/8); % lsr of the family computed as rho(s.l.p.)^(1 / degree(s.l.p.))

A = A/rho; % We consider the normalized family for simplicity and to avoid large entries appearing when we let the algorithms investigate very long products.


% Setting parameters for the algoritms: delta (accuracy), numeric_error = 0 (no numerical error in the entries of A) and display = 1 (to show the progress after each step)
delta = 1e-6; numeric_error = [0 0]; display = 1;



%% Algorithm S (Algorithm 4.1 in the main paper)

M=[500 5000 10000]; % Different values of maximum number of allowed antinorm evaluations for comparison

V_in = eye(2); % 1-antinorm as the polytope antinorm

for i = 1:length(M)
    [lsr_standard(i,:),p_standard(i,:),~] = lsr_comp_standard(A, numeric_error, delta, M(i), V_in, [], display);
end



%% Algorithm A with initial antinorm the 1-antinorm

M = [100 500 1000];

V_adaptive_A = []; % To store the adaptive antinorms obtained at the end of the algorithm for all values of M

for i = 1:length(M)

    [lsr_adaptive_A(i,:),p_adaptive_A(i,:),~,V_new] = adaptive_subradius_comp(A, numeric_error, delta, M(i), V_in, [], display);
    
    V_adaptive_A = [V_new zeros(2,1)]; % The column (0, 0)^T is added to separate the different vertices set obtained for different values of M.

end



%% Algorithm E with initial antinorm the 1-antinorm

M = [50 500 1000];

V_adaptive_E = []; % To store the adaptive antinorms obtained at the end of the algorithm for all values of M

theta = 1.005; % The scaling parameter (to rescale eigenvectors prior to adding them)

for i = 1:length(M)

    [lsr_adaptive_E(i,:),p_adaptive_E(i,:),~,V_new] = adaptive_eigenvectors_subradius_comp(A, numeric_error, delta, M(i), V_in, [], display,theta);
    
    V_adaptive_E = [V_new zeros(2,1)];

end