function [lower,upper,x] = real_antinorm(V,z)

% This function (Algorithm SM.1.1 in the supplementary material) computes the polytope antinorm (corresponding to the vertex set V) of a vector z.


%% Input
% V is a dxp real matrix which contains the p vertices defining the polytope antinorm as columns
% z is a real vector which belongs to the cone on which the antinorm is defined

%% Output
% The values lower and upper provide a lower and upper bound for the value of the antinorm of z
% The vector x consists of 1/a(z) as its first entry, and the values of the variables in the constraints of the LP problem in all other entries



[d,p]=size(V); % d = dimension, p = number of vertices
A=[-z V;0 -ones(1,m)]; % augmented matrix as defined in Algorithm SM.1.1

% The vector c = (1,0,...,0) is used since linprog solves min c^T x, which would be exactly x(1) in this case
c = zeros(m+1,1); c(1,1) = 1;

% The vector b is used to take care of the inequality constraints in the LP problem
b = [zeros(length(z),1);-1];

% Aeq and beq are empty because there are no equality constraints in this LP problem
Aeq=[]; beq=[];

% Tolerance to solve with linprog and other options
tol=1e-10;
options=optimset('MaxIter',max(300,n),'DISPLAY','off','TolFun',tol);

% x solution of the LP problem using linprog, EXITFLAG is also returned as it gives information in case an error should appear
[x,~,EXITFLAG,~,~]=linprog(c,A,b,Aeq,beq,zeros(m+1,1),[],options);

if (EXITFLAG<0) % If EXITFLAG is less than 0, then there is an error
    
    if (EXITFLAG==-2) % when EXITFLAG = -2, no feasible point was found. In this case, we set the antinorm of z equal to 0.
        lower = 0; upper = 0;
    
    else % otherwise (i.e., EXITFLAG = -3,...,-9), we print a message which gives the exact EXITFLAG, allowing us to handle the issue properly 
        disp('Problems with the minimizer. Exit flag: '),disp(EXITFLAG);
        return
    end
    
end

if isempty(x) % if no solution was found, return 0
    lower=0; upper=0;
    return
end

% Computation of a(z) = 1/x(1). For simplicity, when x(1) is +inf, we set a(z)=0
if isinf(x(1))
    lower=0; upper=0;
else
    upper = 1/(x(1)-tol); lower = 1/(x(1)+tol);
end
