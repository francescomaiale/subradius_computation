function [Radius, m, errOut, V, X_opt] = adaptive_subradius_comp(...
    A, dA, delta, M, V, errIn, Display, pruneFlag)
%==========================================================================
% Description:
%   Estimate the spectral subradius (lower spectral radius) of a matrix
%   family F = {A_1, …, A_N} via an adaptive subradius computation
%   using polytope antinorms (Algorithm A).
%
% Inputs (positional or name-value):
%   A         : matrix of size (n × n*N) or (n*N × n), representing
%               N blocks A_i ∈ ℝ^(n×n)
%   dA        : scalar or [N×1] vector of additive error bounds per
%               block (default: eps)
%   delta     : [tol_low; tol_high] tolerances for convergence
%               (default: [1e-10; 1e-10])
%   M         : maximum total products to explore (default: 500)
%   V         : initial vertices defining polytope antinorm
%               (default: eye(n))
%   errIn     : [err_low; err_high] relative error bounds for aNorm
%               (default: [2^-50*n; 2^-49*n])
%   Display   : logical flag to print iteration summaries
%               (default: false)
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
narginchk(1,8);
[ma, na] = size(A);

if nargin < 2 || isempty(dA),        dA = [];            end
if nargin < 3 || isempty(delta),     delta = [1e-10;1e-10]; end
if nargin < 4 || isempty(M),         M     = 500;        end
if nargin < 5 || isempty(V),         V     = eye(ma);    end
if nargin < 6 || isempty(errIn),     errIn = [];         end
if nargin < 7 || isempty(Display),   Display = false;    end
if nargin < 8 || isempty(pruneFlag), pruneFlag = true;   end

delta = delta(:);
if numel(delta)==1, delta(2)=delta(1); end
assert(isnumeric(delta) && numel(delta)==2, ...
    'delta must be a scalar or 2-element vector.');

tol          = 1e-10;
tol_vertices = 1e-10;

%% 2. Block‐structure check
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

%% 3. Error‐bounds defaults & validation
if isempty(errIn)
    errIn = [2^(-50)*n; 2^(-49)*n];
end
assert(isnumeric(errIn) && numel(errIn)==2, ...
    'errIn must be a 2-element vector.');

%% 4. Prepare per‐block additive errors
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
LowRadius    = 0;
NORM         = zeros(N,1);
LowBound     = zeros(N,1);
NormErr      = zeros(N,1);
v_cand       = zeros(ma,N);
cand_counter = 1;
bestUpperVal = Inf;
X_opt        = [];

%% 6. First‐level products
for i = 1:N
    Y = A(:,(i-1)*n+1:i*n);
    [aN,v_new]   = aNorm(Y,V);
    NORM(i)      = aN*(1+errIn(2));
    LowBound(i)  = NORM(i);
    if proportional
        NormErr(i) = dA(1)*NORM(i);
    else
        NormErr(i) = dA(i);
    end
    rho_Y        = abs(eigs(Y,1,'largestabs'));
    HighRadius   = min(HighRadius, rho_Y);
    if rho_Y < bestUpperVal
        bestUpperVal = rho_Y;
        X_opt        = Y;
    end
    X = [X; Y];
    [~,upper,~] = real_antinorm(V,v_new);
    if upper<=1+tol && norm(v_new)>1e-3
        v_cand(:,cand_counter) = v_new;
        cand_counter = cand_counter+1;
    end
end

v_cand = v_cand(:,1:cand_counter-1);
V = [V, v_cand];
if pruneFlag, V = pruneVertices(V,tol_vertices); end

dA        = NormErr;
NormA     = NORM;
LowRadius = min(NORM + NormErr);
mm        = N;
mCount    = 1;
JJ        = N;
MaxJJ     = N;
BestPower = 1;
ell_opt   = 1;

%% 7. Main adaptive loop
while mm < M && LowRadius < HighRadius - delta(2)
    mCount     = mCount + 1;
    cand_counter = 0;
    v_cand     = zeros(ma,0);
    
    OldHighRadius = HighRadius;
    OldLowRadius  = LowRadius;
    LowRadius     = HighRadius;
    NewJJ         = 0;
    
    newCap       = JJ * N;
    NewNorm      = zeros(newCap,1);
    NewNormErr   = zeros(newCap,1);
    NewLowBound  = zeros(newCap,1);
    NewX         = zeros(n*newCap, n);
    newIndex     = 0;
    
    for k = 1:JJ
        Xk = X((k-1)*n+1:k*n, :);
        for i = 1:N
            newIndex = newIndex + 1;
            
            XX = Xk * A(:, (i-1)*n+1:i*n);
            
            % propagate errors
            NewNormErr(newIndex) = ...
                NORM(k)*(errIn(1)*NormA(i) + dA(i)) + ...
                NormErr(k)*(NormA(i)       + dA(i));
            
            [aN, v_new]            = aNorm(XX, V);
            NewNorm(newIndex)      = aN*(1+errIn(2));
            
            [~,upper,~] = real_antinorm(V,v_new);
            if upper<=1+tol && norm(v_new)>1e-3
                cand_counter = cand_counter+1;
                v_cand(:,cand_counter) = v_new;
            end
            
            candidate = (NewNorm(newIndex) + NewNormErr(newIndex))^(1/mCount);
            NewLowBound(newIndex) = max(LowBound(k), candidate);
            
            rho_XX         = abs(eigs(XX,1,'largestabs'));
            currentUpperVal= rho_XX^(1/mCount);
            if currentUpperVal < bestUpperVal
                bestUpperVal = currentUpperVal;
                X_opt        = XX;
            end
            HighRadius = min(HighRadius, currentUpperVal);
            
            if NewLowBound(newIndex) < OldHighRadius - delta(1)
                NewJJ = NewJJ + 1;
                LowRadius = min(LowRadius, NewLowBound(newIndex));
                NewX((NewJJ-1)*n+1:NewJJ*n, :) = XX;
                % Compact metadata with its retained product. Discarded
                % products must not supply bounds to the next level.
                NewNorm(NewJJ) = NewNorm(newIndex);
                NewNormErr(NewJJ) = NewNormErr(newIndex);
                NewLowBound(NewJJ) = NewLowBound(newIndex);
            end
        end
    end
    
    % update counters and queues
    mm      = mm + JJ*N;
    NORM    = NewNorm(1:NewJJ);
    LowBound= NewLowBound(1:NewJJ);
    NormErr = NewNormErr(1:NewJJ);
    JJ      = NewJJ;
    MaxJJ   = max(JJ, MaxJJ);

    % Account for the tolerance used to discard products (Algorithm 6,
    % line 31), including when the retained queue is empty.
    LowRadius = max(OldLowRadius, min(LowRadius, HighRadius - delta(1)));
    
    if (HighRadius - LowRadius) < (OldHighRadius - OldLowRadius)
        BestPower = mCount;
    end
    if HighRadius < OldHighRadius
        ell_opt = mCount;
    end
    
    if Display
        disp([BestPower, mCount, mm, JJ, ell_opt]);
        disp([LowRadius, HighRadius]);
    end
    
    % update vertices
    V = [V, v_cand];
    if pruneFlag, V = pruneVertices(V, tol_vertices); end
    
    if NewJJ == 0
        disp('Terminating early: no new products to explore.');
        break;
    else
        X = NewX(1:NewJJ*n, :);
    end
end

%% 8. Final outputs
[~, idx] = unique(V','rows','stable');
V        = V(:, idx);
Radius   = [LowRadius, HighRadius];
m        = [BestPower, mCount, mm, MaxJJ, ell_opt];
errOut   = errIn;

end


%% Helper: pruneVertices
function Vp = pruneVertices(V, tol)
% Remove vertices v for which the lower antinorm ≥ 1+tol when omitted
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
