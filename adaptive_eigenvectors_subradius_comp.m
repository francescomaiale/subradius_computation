function [Radius,m,Err,V,X_opt] = adaptive_eigenvectors_subradius_comp( ...
    A, dA, delta, M, V, Err, Display, theta, pruneFlag)

if nargin < 2 || isempty(dA),       dA = [];          end
if nargin < 3 || isempty(delta),    delta = 1e-6;     end
if nargin < 4 || isempty(M),        M = 500;          end
if nargin < 5 || isempty(V),        V = eye(size(A,1)); end
if nargin < 6 || isempty(Err),      Err = [];         end
if nargin < 7 || isempty(Display),  Display = 0;      end
if nargin < 8 || theta < 1,         theta = 1.01;     end
if nargin < 9 || isempty(pruneFlag), pruneFlag = true; end

delta = delta(:);
if numel(delta) < 2
    delta = [delta(1); delta(1)];
end

[ma, na] = size(A);
tol_vertices = 1e-8;  % Tolerance for vertex pruning
tol = 1e-8;

if (na > ma && rem(na,ma) == 0)
    N = na / ma;
    n = ma;
elseif (ma > na && rem(ma,na) == 0)
    A = A.';
    [ma, na] = size(A);
    N = na / ma;
    n = ma;
else
    disp('The matrix A has incompatible dimensions.');
    Radius = [0,0];
    m      = [];
    X_opt  = [];
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

bestUpperVal = Inf;
X_opt        = [];

NORM     = zeros(N,1);
LowBound = zeros(N,1);
NormErr  = zeros(N,1);

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

    try
        [~, upperVal, ~] = real_antinorm(V, v_new);
        if upperVal <= 1 + tol
            if norm(v_new) > tol
                V = [V, v_new];
            end
        end
    catch ME
        fprintf('Error on block %d: %s\n', i, ME.message);
        continue;
    end

    rho_Y = abs(eigs(Y, 1, 'largestabs'));

    if HighRadius < tol
        HighRadius = rho_Y;
    else
        HighRadius = min(HighRadius, rho_Y);
    end

    if abs(HighRadius - rho_Y) < tol
        [vY, ~] = eigs(Y, 1, 'largestabs');
        vY = abs(vY);
        if norm(vY) > tol
            [lower_ev,~,~] = real_antinorm(V, vY);
            V = [V, vY/(lower_ev * theta)];
        end
    end

    currentUpperVal = rho_Y;
    if currentUpperVal < bestUpperVal
        bestUpperVal = currentUpperVal;
        X_opt        = Y;
    end

    X = [X; Y];
end

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

%% Main loop
while (mm < M) && (LowRadius < HighRadius - delta(2))

    mCount = mCount + 1;
    OldHighRadius = HighRadius;
    OldLowRadius  = LowRadius;

    LowRadius = HighRadius;
    NewJJ     = 0;

    newCap      = JJ*N;
    NewNorm     = zeros(newCap,1);
    NewNormErr  = zeros(newCap,1);
    NewLowBound = zeros(newCap,1);
    NewX        = zeros(n*newCap, n);
    newIndex    = 0;

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

            [~, upperVal, ~] = real_antinorm(V, v_new);
            if upperVal <= 1 + tol
                if norm(v_new) > tol
                    V = [V, v_new];
                end
            end

            cand = (NewNorm(newIndex) + NewNormErr(newIndex))^(1/mCount);
            NewLowBound(newIndex) = max(LowBound(k), cand);

            rho_XX = abs(eigs(XX, 1, 'largestabs'));
            currHR = rho_XX^(1/mCount);
            HighRadius = min(HighRadius, currHR);

            if abs(HighRadius - currHR) < tol
                [v, ~] = eigs(XX, 1, 'largestabs');
                v = abs(v);
                [lowEv,~,~] = real_antinorm(V, v);
                V = [V, v/(lowEv * theta)];
            end

            if currHR < bestUpperVal
                bestUpperVal = currHR;
                X_opt        = XX;
            end

            if NewLowBound(newIndex) < HighRadius - delta(1)
                NewJJ = NewJJ + 1;
                LowRadius = min(LowRadius, NewLowBound(newIndex));
                NewX((NewJJ-1)*n+1 : NewJJ*n, :) = XX;
            end
        end
    end

    X        = NewX(1 : NewJJ*n, :);
    NORM     = NewNorm(1 : NewJJ);
    LowBound = NewLowBound(1 : NewJJ);
    NormErr  = NewNormErr(1 : NewJJ);
    JJ       = NewJJ;
    MaxJJ    = max(JJ, MaxJJ);

    LowRadius = max(OldLowRadius, min(LowRadius, HighRadius - delta(1)));
    mm        = mm + JJ*N;

    if NewJJ > 0
        if (HighRadius - LowRadius) < (OldHighRadius - OldLowRadius)
            BestPower = mCount;
        end
    end

    if HighRadius < OldHighRadius
        ell_opt = mCount;
    end

    if Display
        disp([BestPower, mCount, mm, JJ, ell_opt]);
        disp([LowRadius, HighRadius]);
    end

    if pruneFlag
        V = pruneVertices(V, tol_vertices);
    end

    if NewJJ == 0
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


function Vp = pruneVertices(V, tol_vertices)

if norm(V) < 1e-3
    Vp = V;
else

    reduce = false;
    while ~reduce
        n_v = size(V,2);
        if n_v <= 1
            break;
        end
        reduce = true;
        for i = n_v : -1 : 1
            v = V(:, i);
            W = [V(:,1:i-1), V(:,i+1:n_v)];
            if rank(W) == 1
                break;
            end
            [lowerVal, ~, ~] = real_antinorm(W, v);
            if lowerVal >= 1 + tol_vertices
                reduce = false;
                V = W;
                break;
            end
        end
    end
    [~, idx] = unique(V', 'rows', 'stable');
    Vp = V(:, idx);
end
end
