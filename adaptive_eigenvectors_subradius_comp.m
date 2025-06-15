function [Radius, m, errOut, V, X_opt] = adaptive_eigenvectors_subradius_comp(...
    A, dA, delta, M, V, errIn, Display, theta, pruneFlag)
%==========================================================================
% Description:
%   Estimate the spectral subradius (lower spectral radius) of a matrix
%   family F = {A_1, …, A_N} using an adaptive eigenvector‐based
%   subradius computation with polytope antinorms and eigenvector expansion.
%
% Inputs (positional or name-value):
%   A         : matrix of size (n × n*N) or (n*N × n), representing
%               N blocks A_i ∈ ℝ^(n×n)
%   dA        : scalar or [N×1] vector of additive error bounds per
%               block (default: eps)
%   delta     : [tol_low; tol_high] tolerances for convergence
%               (default: [1e-6; 1e-6])
%   M         : maximum total products to explore (default: 500)
%   V         : initial vertices defining polytope antinorm
%               (default: eye(n))
%   errIn     : [err_low; err_high] relative error bounds for aNorm
%               (default: [2^-50*n; 2^-49*n])
%   Display   : logical flag to print iteration summaries
%               (default: false)
%   theta     : expansion factor for eigenvector‐based vertices
%               (default: 1.01)
%   pruneFlag : logical flag to prune redundant vertices
%               (default: true)
%
% Outputs:
%   Radius    : [LowRadius, HighRadius] bounds on the true subradius
%   m         : [BestPower, IterCount, TotalEvals, MaxQueue, SLPPower]
%   errOut    : [err_low, err_high] actual error bounds used
%   V         : final set of vertices
%   X_opt     : product achieving the current best upper bound
%
% NOTE:
%   Relies on aNorm and real_antinorm. It can be changed to the CVX
%   solver (real_antinorm_cvx) by following the instructions in aNorm.m 
%   Pruning uses pruneVertices below.
%==========================================================================

%% 1. Input parsing & defaults
narginchk(1,9);
[ma, na] = size(A);

if nargin < 2 || isempty(dA),        dA = [];             end
if nargin < 3 || isempty(delta),     delta = [1e-6;1e-6]; end
if nargin < 4 || isempty(M),         M     = 500;         end
if nargin < 5 || isempty(V),         V     = eye(ma);     end
if nargin < 6 || isempty(errIn),     errIn = [];          end
if nargin < 7 || isempty(Display),   Display = false;     end
if nargin < 8 || isempty(theta) || theta<1, theta = 1.01;  end
if nargin < 9 || isempty(pruneFlag), pruneFlag = true;    end

delta = delta(:);
if numel(delta)==1, delta(2)=delta(1); end
assert(isnumeric(delta) && numel(delta)==2, ...
    'delta must be a scalar or 2-element vector.');

tol          = 1e-8;
tol_vertices = 1e-8;

%% 2. Block-structure check
if na>ma && rem(na,ma)==0
    N = na/ma;  n = ma;
elseif ma>na && rem(ma,na)==0
    A = A.';  [ma, na] = size(A);
    N = na/ma;  n = ma;
else
    disp('The matrix A has incompatible dimensions. Aborting...');
    Radius = [0,0];  m = [];  errOut = [];  V = [];  X_opt = [];
    return;
end

%% 3. Error-bounds defaults & validation
if isempty(errIn)
    errIn = [2^(-50)*n; 2^(-49)*n];
end
assert(isnumeric(errIn) && numel(errIn)==2, ...
    'errIn must be a 2-element vector.');

%% 4. Prepare per-block additive errors
dA = dA(:);
if numel(dA)==n
    proportional = false;
else
    proportional = true;
    if isempty(dA), dA = eps; end
end

%% 5. Initialization
X            = [];
HighRadius   = Inf;
bestUpperVal = Inf;
X_opt        = [];

NORM         = zeros(N,1);
LowBound     = zeros(N,1);
NormErr      = zeros(N,1);

%% 6. First-level products
for i = 1:N
    Y = A(:, (i-1)*n+1 : i*n);

    [aN, v_new]   = aNorm(Y, V);
    NORM(i)       = aN*(1 + errIn(2));
    LowBound(i)   = NORM(i);

    if proportional
        NormErr(i) = dA(1)*NORM(i);
    else
        NormErr(i) = dA(i);
    end

    % include new vertex if antinorm upper ≤ 1
    try
        [~, upperVal, ~] = real_antinorm(V, v_new);
        if upperVal <= 1 + tol && norm(v_new)>tol
            V = [V, v_new];
        end
    catch ME
        fprintf('Warning: real_antinorm failed on block %d: %s\n', i, ME.message);
    end

    % spectral radius and eigenvector expansion
    rho_Y = abs(eigs(Y,1,'largestabs'));
    HighRadius = min(HighRadius, rho_Y);

    if abs(HighRadius - rho_Y) < tol
        [vY, ~] = eigs(Y,1,'largestabs');
        vY = abs(vY);
        [lower_ev,~,~] = real_antinorm(V, vY);
        V = [V, vY/(lower_ev*theta)];
    end

    if rho_Y < bestUpperVal
        bestUpperVal = rho_Y;
        X_opt        = Y;
    end

    X = [X; Y];
end

%% prepare loop variables
dA        = NormErr;
NormA     = NORM;
LowRadius = min(NORM + NormErr);
mm        = N;
mCount    = 1;
JJ        = N;
MaxJJ     = N;
BestPower = 1;
ell_opt   = 1;

if pruneFlag
    V = pruneVertices(V, tol_vertices);
end

%% 7. Main adaptive loop
while mm < M && LowRadius < HighRadius - delta(2)
    mCount      = mCount + 1;
    OldHigh     = HighRadius;
    OldLow      = LowRadius;
    LowRadius   = HighRadius;
    NewJJ       = 0;

    cap         = JJ * N;
    NewNorm     = zeros(cap,1);
    NewNormErr  = zeros(cap,1);
    NewLowBound = zeros(cap,1);
    NewX        = zeros(n*cap, n);
    idx         = 0;

    for k = 1:JJ
        Xk = X((k-1)*n+1:k*n, :);
        for i = 1:N
            idx = idx + 1;
            XX  = Xk * A(:, (i-1)*n+1:i*n);

            % propagate error
            NewNormErr(idx) = ...
                NORM(k)*(errIn(1)*NormA(i) + dA(i)) + ...
                NormErr(k)*(NormA(i)        + dA(i));

            [aN, v_new]      = aNorm(XX, V);
            NewNorm(idx)     = aN*(1 + errIn(2));

            % antinorm upper test
            [~, upperVal, ~] = real_antinorm(V, v_new);
            if upperVal <= 1 + tol && norm(v_new)>tol
                V = [V, v_new];
            end

            candVal = (NewNorm(idx) + NewNormErr(idx))^(1/mCount);
            NewLowBound(idx) = max(LowBound(k), candVal);

            rho_XX   = abs(eigs(XX,1,'largestabs'))^(1/mCount);
            HighRadius = min(HighRadius, rho_XX);

            % spectral‐expansion test
            if abs(HighRadius - rho_XX) < tol
                [v,~] = eigs(XX,1,'largestabs');
                v = abs(v);
                [lowEv,~,~] = real_antinorm(V, v);
                V = [V, v/(lowEv*theta)];
            end

            if rho_XX < bestUpperVal
                bestUpperVal = rho_XX;
                X_opt        = XX;
            end

            if NewLowBound(idx) < OldHigh - delta(1)
                NewJJ = NewJJ + 1;
                LowRadius = min(LowRadius, NewLowBound(idx));
                NewX((NewJJ-1)*n+1:NewJJ*n, :) = XX;
            end
        end
    end

    % update queues & counters
    mm        = mm + JJ*N;
    NORM      = NewNorm(1:NewJJ);
    LowBound  = NewLowBound(1:NewJJ);
    NormErr   = NewNormErr(1:NewJJ);
    JJ        = NewJJ;
    MaxJJ     = max(MaxJJ, JJ);

    % update best power and SLP power
    if (HighRadius - LowRadius) < (OldHigh - OldLow)
        BestPower = mCount;
    end
    if HighRadius < OldHigh
        ell_opt = mCount;
    end

    if Display
        disp([BestPower, mCount, mm, JJ, ell_opt]);
        disp([LowRadius, HighRadius]);
    end

    if pruneFlag
        V = pruneVertices(V, tol_vertices);
    end

    if JJ == 0
        disp('Terminating early: no new products to explore.');
        break;
    else
        X = NewX(1:JJ*n, :);
    end
end

%% 8. Final outputs
[~, idxv] = unique(V','rows','stable');
V         = V(:,idxv);
Radius    = [LowRadius, HighRadius];
m         = [BestPower, mCount, mm, MaxJJ, ell_opt];
errOut    = errIn;

end

%% Helper: pruneVertices
function Vp = pruneVertices(V, tol)
% Remove vertices v for which lower antinorm ≥ 1+tol when omitted
    changed = true;
    while changed && size(V,2)>1
        changed = false;
        for i = size(V,2):-1:1
            v = V(:,i);
            W = [V(:,1:i-1), V(:,i+1:end)];
            if rank(W)==1, continue; end
            [lowerVal,~,~] = real_antinorm(W, v);
            if lowerVal >= 1+tol
                V = W;
                changed = true;
                break;
            end
        end
    end
    [~, idx] = unique(V','rows','stable');
    Vp       = V(:,idx);
end
