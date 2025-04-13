function [Radius, m, Err, V, X_opt] = adaptive_subradius_comp(A, dA, delta, M, V, Err, Display, pruneFlag)

if nargin < 2 || isempty(dA),    dA = [];        end
if nargin < 3 || isempty(delta), delta = 1e-10;   end
if nargin < 4 || isempty(M),     M = 500;        end
if nargin < 5 || isempty(V),     V = eye(size(A,1)); end
if nargin < 6 || isempty(Err),   Err = [];       end
if nargin < 7 || isempty(Display), Display = 0;  end
if nargin < 8 || isempty(pruneFlag), pruneFlag = true; end

delta = delta(:);
if numel(delta) < 2
    delta = [delta(1); delta(1)];
end

[ma, na] = size(A);
tol = 1e-10;
tol_vertices = 1e-10;

if na > ma && rem(na, ma) == 0
    N = na / ma;
    n = ma;
elseif ma > na && rem(ma, na) == 0
    A = A.';
    [ma, na] = size(A);
    N = na / ma;
    n = ma;
else
    disp('The matrix A has incompatible dimensions.');
    Radius = [0, 0];
    m = [];
    X_opt = [];
    return;
end

if numel(Err) < 2
    Err = [2^(-50)*n, 2^(-49)*n];
end

dA = dA(:);
if numel(dA) == n
    proportional = 0;
else
    proportional = 1;
    if isempty(dA)
        dA = eps;
    end
end

%% First Step

X = [];
HighRadius = Inf;
LowRadius  = 0;

NORM     = zeros(N,1);
LowBound = zeros(N,1);
NormErr  = zeros(N,1);
v_cand = zeros(ma,N);
cand_counter = 1;

bestUpperVal = Inf;
X_opt = [];

for i = 1 : N
    Y = A(:, (i-1)*n + 1 : i*n);

    [aN, v_new] = aNorm(Y, V);
    NORM(i)     = aN * (1 + Err(2));
    LowBound(i) = NORM(i);

    if proportional
        NormErr(i) = dA(1) * NORM(i);
    else
        NormErr(i) = dA(i);
    end

    rho_Y = abs(eigs(Y, 1, 'largestabs'));
    HighRadius = min(HighRadius, rho_Y);

    if rho_Y < bestUpperVal
        bestUpperVal = rho_Y;
        X_opt = Y;
    end

    X = [X; Y];

    [~, upper, ~] = real_antinorm(V, v_new);
    if upper <= 1 + tol
        %V = [V, v_new];
        if norm(v_new) > 1e-3
        v_cand(:,cand_counter) = v_new;
        cand_counter = cand_counter + 1;
        end
    end
end

v_cand = v_cand(:,1:cand_counter-1);
V = [V v_cand];

if pruneFlag
    V = pruneVertices(V, tol_vertices);
end

dA       = NormErr;
NormA    = NORM;
LowRadius= min(NORM + NormErr);
mm       = N;
mCount   = 1;
JJ       = N;
MaxJJ    = N;
BestPower= 1;
ell_opt  = 1;

%% Main loop
while (mm < M) && (LowRadius < HighRadius - delta(2))

    mCount = mCount + 1;
    cand_counter = 0;
    v_cand = [];

    OldHighRadius = HighRadius;
    OldLowRadius  = LowRadius;
    LowRadius     = HighRadius;
    NewJJ         = 0;

    newCap       = JJ*N;
    NewNorm      = zeros(newCap, 1);
    NewNormErr   = zeros(newCap, 1);
    NewLowBound  = zeros(newCap, 1);
    NewX         = zeros(n*newCap, n);
    newIndex     = 0;

    for k = 1 : JJ
        Xk = X((k-1)*n+1 : k*n, :);
        for i = 1 : N
            newIndex = newIndex + 1;

            XX = Xk * A(:, (i-1)*n + 1 : i*n);

            NewNormErr(newIndex) = ...
                ( NORM(k)    * (Err(1)*NormA(i) + dA(i)) + ...
                NormErr(k) * (NormA(i)        + dA(i)) );

            [aN, v_new] = aNorm(XX, V);
            NewNorm(newIndex) = aN * (1 + Err(2));

            [~, upper, ~] = real_antinorm(V, v_new);
            if upper <= 1 + tol
                %V = [V, v_new];
                if norm(v_new) > 1e-3
                cand_counter = cand_counter+1;
                v_cand = [v_cand v_new];
                end
            end

            candidate = (NewNorm(newIndex) + NewNormErr(newIndex))^(1/mCount);
            NewLowBound(newIndex) = max(LowBound(k), candidate);

            rho_XX = abs(eigs(XX, 1, 'largestabs'));
            currentUpperVal = rho_XX^(1/mCount);
            if currentUpperVal < bestUpperVal
                bestUpperVal = currentUpperVal;
                X_opt        = XX;
            end
            HighRadius = min(HighRadius, currentUpperVal);

            if NewLowBound(newIndex) < OldHighRadius - delta(1)
                NewJJ = NewJJ + 1;
                LowRadius = min(LowRadius, NewLowBound(newIndex));
                NewX((NewJJ-1)*n + 1 : NewJJ*n, :) = XX;
            end

        end
    end

    LowRadius = max(OldLowRadius, min(LowRadius, HighRadius - delta(1)));

    mm      = mm + JJ*N;
    NORM    = NewNorm(1 : NewJJ);
    LowBound= NewLowBound(1 : NewJJ);
    NormErr = NewNormErr(1 : NewJJ);
    JJ      = NewJJ;
    MaxJJ   = max(JJ, MaxJJ);

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

    V = [V v_cand];

    if pruneFlag
        V = pruneVertices(V, tol_vertices);
    end

    if NewJJ == 0
        disp('int');
        break;
    else
        X = NewX(1 : NewJJ*n, :);
    end

end

[~, idx] = unique(V', 'rows', 'stable');
V = V(:, idx);

m      = [BestPower, mCount, mm, MaxJJ, ell_opt];
Radius = [LowRadius, HighRadius];
end

function Vp = pruneVertices(V, tol)
reduce = false;
while ~reduce
    n_v = size(V,2);
    if n_v <= 1
        break;
    end
    reduce = true;
    for i = n_v : -1 : 1
        v = V(:, i);
        W = [V(:, 1:i-1), V(:, i+1:n_v)];

        if rank(W) == 1
            break
        end

        [lowerVal, ~, ~] = real_antinorm(W, v);
        if lowerVal >= 1 + tol
            reduce = false;
            V = W;
            break;
        end
    end
end
[~, idx] = unique(V', 'rows', 'stable');
Vp = V(:, idx);
end
