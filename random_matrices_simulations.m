% In this file, the strategy proposed in Section 6.2 of the main paper is applied to randomly generated families of full and sparse matrices (with both adaptive algorithms).
% Depending on the matrix dimensions, the maximum number of iterations (MaxIter), the precision delta, and other parameters (e.g., M), the simulations may take a considerable amount of time to complete.
% To get a better understanding of the computational time required for subradius computation, it is recommended to start by commenting out either the full or sparse matrix sections and setting all parameters for a single value of the dimension d.
% This approach will provide a clearer idea of the time-consuming nature of the subradius computation before running the entire set of simulations.


% This script runs simulations for all values of the dimension d selected below. It is recommended to keep d <= 200, particularly when working with sparse matrices.
%dimension_values = [25 50 100 150 200];
dimension_values = 50;


%% Simulations for randomly generated full matrices

% Preallocate space for lower and upper bounds, performance metric, computational time, and number of vertices of the adaptive antinorm
lower_bound_full = zeros(length(dimension_values),1);
upper_bound_full = zeros(length(dimension_values),1);
time_full = zeros(length(dimension_values),1);
performance_metric_full = zeros(length(dimension_values),5);
vertices_number_full = zeros(length(dimension_values),1);
Iter = zeros(length(dimension_values),1); % tracks how many iterations the procedure does; if Iter = MaxIter, convergence to delta has likely not been achieved

for ell = 1:length(dimension_values)
    
    tic % start the timer
    d = dimension_values(ell); % current dimension
    A_1 = rand(d,d); A_2 = rand(d,d); A = [A_1 A_2];
    A_in = A; % store the initial family since it will be rescaled multiple times during the procedure
    
   numeric_error = [0 0]; % assume that the entries of A_1 and A_2 are known with no numerical error
   
    % As the initial antinorm, we take the one given by the leading eigenvector of A_1
    [v,d] = eigs(A_1,1,'largestabs'); v = abs(v(:));
    V_in = v;
    
    % Apply either Algorithm A or E with a small value of M to obtain a preliminary bound and use the lower bound to rescale the family.
    M = 10; theta = 1.005; delta = 1e-6;
    [preliminary_bounds,~,~,~] = adaptive_subradius_comp(A,numeric_error,delta,M,V_in,[],0);
    
    % Update the lower and upper bounds with the preliminary bound
    lower_bound_full(ell) = preliminary_bounds(1);
    upper_bound_full(ell) = preliminary_bounds(2);

    A = A_in / preliminary_bounds(1); % rescaling of the family
    M = 500; % this is the value of M that is used for all other iterations
    display = 0; % does not show progression
    MaxIter = 10; % if convergence is not achieved, this terminates the algorithm after a certain number of iterations
    
    CurrentIter = 0; % current number of iterations (0)
    
    % We introduce lsr(1) = lower bound and lsr(2) = upper bound because we want to use them to stop the iterations once lsr(2)-lsr(1) is small enough
    lsr(2)=preliminary_bounds(2) / preliminary_bounds(1); 
    lsr(1)=1;
    
    V = V_in; % V is used to store the adaptive antinorm, since each iteration will start from the antinorm refined in the previous one

    while (lsr(2) - lsr(1) >= delta && CurrentIter <= MaxIter)
    
        CurrentIter = CurrentIter + 1;
        disp('Start iteration '); disp(CurrentIter);

        [lsr,perf_metric,~,V] = adaptive_eigenvectors_subradius_comp(A,numeric_error,delta,M,V,[],display,theta);
        % [lsr,perf_metric,~,V] = adaptive_subradius_comp(A,numeric_error,delta,M,V,[],display);
        
        lower_bound_full(ell) = lower_bound_full(ell) * lsr(1); % update the lower bound
        upper_bound_full(ell) = lower_bound_full(ell) * lsr(2) / lsr(1); % update the upper bound
        A = A / lsr(1); % rescale the family further
        disp('End iteration '); disp(CurrentIter);
    end

    performance_metric_full(ell,:) = perf_metric(:);
    Iter(ell) = CurrentIter(ell);
    vertices_number_full(ell)=size(V,2);
    time_full(ell)=toc; % stop the timer and save the time
end



%% Simulations for randomly generated sparse matrices

sparsity_density = 0.25; % select the sparsity density between 0.2 and 1.0; to explore lower densities (e.g., 0.1), perturbation theory must be applied as in Section 3.3 in the main paper

% Preallocate space for lower and upper bounds, performance metric, computational time, and number of vertices of the adaptive antinorm
lower_bound_sparse = zeros(length(dimension_values),1);
upper_bound_sparse = zeros(length(dimension_values),1);
time_sparse = zeros(length(dimension_values),1);
performance_metric_sparse = zeros(length(dimension_values),5);
vertices_number_sparse = zeros(length(dimension_values),1);
Iter_sparse = zeros(length(dimension_values),1); % tracks how many iterations the procedure does; if Iter = MaxIter, convergence to delta has likely not been achieved

for ell = 1:length(dimension_values)
    
    tic % start the timer
    d = dimension_values(ell); % current dimension
    A_1 = sprand(d,d,sparsity_density); A_2 = sprand(d,d,sparsity_density); 
    A = [A_1 A_2];
    A_in = A; % store the initial family since it will be rescaled multiple times during the procedure
    
   numeric_error = [0 0]; % assume that the entries of A_1 and A_2 are known with no numerical error
   
    % As the initial antinorm, we take the one given by the leading eigenvector of A_1
    [v,d] = eigs(A_1,1,'largestabs'); v = abs(v(:));
    V_in = v;
    
    % Apply either Algorithm A or E with a small value of M to obtain a preliminary bound and use the lower bound to rescale the family.
    M = 10; theta = 1.005; delta = 1e-6;
    [preliminary_bounds,~,~,~] = adaptive_subradius_comp(A,numeric_error,delta,M,V_in,[],0);
    
    % Update the lower and upper bounds with the preliminary bound
    lower_bound_sparse(ell) = preliminary_bounds(1);
    upper_bound_sparse(ell) = preliminary_bounds(2);

    A = A_in / preliminary_bounds(1); % rescaling of the family
    M = 500; % this is the value of M that is used for all other iterations
    display = 0; % does not show progression
    MaxIter = 10; % if convergence is not achieved, this terminates the algorithm after a certain number of iterations
    
    CurrentIter = 0; % current number of iterations (0)
    
    % We introduce lsr(1) = lower bound and lsr(2) = upper bound because we want to use them to stop the iterations once lsr(2)-lsr(1) is small enough
    lsr(2) = preliminary_bounds(2) / preliminary_bounds(1); 
    lsr(1) = 1;
    
    V = V_in; % V is used to store the adaptive antinorm, since each iteration will start from the antinorm refined in the previous one

    while (lsr(2) - lsr(1) >= delta && CurrentIter <= MaxIter)
    
        CurrentIter = CurrentIter + 1;
        disp('Start iteration '); disp(CurrentIter);

        [lsr,perf_metric,~,V] = adaptive_eigenvectors_subradius_comp(A,numeric_error,delta,M,V,[],display,theta);
        % [lsr,perf_metric,~,V] = adaptive_subradius_comp(A,numeric_error,delta,M,V,[],display);
        
        lower_bound_sparse(ell) = lower_bound_sparse(ell) * lsr(1); % update the lower bound
        upper_bound_sparse(ell) = lower_bound_sparse(ell) * lsr(2) / lsr(1); % update the upper bound
        A = A / lsr(1); % rescale the family further
        disp('End iteration '); disp(CurrentIter);
    end

    performance_metric_sparse(ell,:) = perf_metric(:);
    Iter_sparse(ell) = CurrentIter(ell);
    vertices_number_sparse(ell)=size(V,2);
    time_sparse(ell)=toc; % stop the timer and save the time
end