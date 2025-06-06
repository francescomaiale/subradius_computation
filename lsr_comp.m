function [lsr, m, errOut] = lsr_comp(A, dA, delta, M, V, errIn, Display)
%==========================================================================
% Description:
%   Estimate the lower spectral radius (LSR) of a matrix family
%   F = {A_1, …, A_N} by the fixed antinorm algorithm (Algorithm S).
%
% Inputs (positional or name-value):
%   A       : matrix of size (n × n*N) or (n*N × n), representing N blocks A_i ∈ ℝ^(n×n)
%   dA      : scalar or [N×1] vector of additive error bounds per block (default: zeros(N,1))
%   delta   : [tol_low; tol_high] tolerances for convergence (default: [1e-6; 1e-6])
%   M       : max total antinorm evaluations (default: 500)
%   V       : vertices defining polytope antinorm (default: eye(n), the 1-antinorm)
%   errIn   : [err_low; err_high] relative error bounds for aNorm (default: [2^-50*n; 2^-49*n])
%   Display : logical flag to print iteration summary (default: false)
%
% Outputs:
%   lsr     : [LowRadius, HighRadius] bounds on the true LSR
%   m       : [BestPower; FinalDegree; TotalEvals; MaxQueue; SLPPower]
%             • BestPower    – degree k minimizing (High_k − Low_k)
%             • FinalDegree  – last degree k explored by the algorithm
%             • TotalEvals   – total antinorm evaluations performed
%             • MaxQueue     – maximum number of products retained at any k
%             • SLPPower     – degree at which High_k strictly decreases
%   errOut  : [err_low, err_high] actual error bounds used
%==========================================================================

% 1. Input parsing & defaults
narginchk(1,7);
[ma, na] = size(A);

if nargin < 2 || isempty(dA),      dA = [];      end
if nargin < 3 || isempty(delta),    delta = [1e-6;1e-6]; end
if nargin < 4 || isempty(M),        M = 500;      end
if nargin < 5 || isempty(V),        V = eye(ma);  end
if nargin < 6 || isempty(errIn),    errIn = [];   end
if nargin < 7 || isempty(Display),  Display = false; end

delta = delta(:);
if numel(delta)==1, delta(2)=delta(1); end
assert(isnumeric(delta) && numel(delta)==2, 'delta must be scalar or 2-element vector.');


    % Attempt to interpret if A is (n x n*N) or (n*N x n)
    if na > ma && rem(na, ma) == 0
        N = na / ma;  % number of sub-blocks
        n = ma;       % each block is (n x n)
    elseif ma > na && rem(ma, na) == 0
        A = A.';
        [ma, na] = size(A);
        N = na / ma;
        n = ma;
    else
        disp('The matrix A has incompatible dimensions. Aborting...');
        lsr = [0, 0]; 
        m   = []; 
        return;
    end

    if isempty(err)
        err = [2^(-50)*n, 2^(-49)*n];
    end
    assert( numel(err) == 2, '''err'' must be a 2-element vector.' );

    dA = dA(:);
    if numel(dA) == n
        prop = 0;  
    else
        prop = 1; 
        if isempty(dA)
            dA = eps;
        end
    end

    X = [];  

    HighRadius = inf;

    normVals  = zeros(N,1); 
    LowBound  = zeros(N,1); 
    NormErr   = zeros(N,1);

    for i = 1:N
        Y = A(:, (i-1)*n + 1 : i*n);

        % Compute antinorm of Y wrt V:
        [aN, ~] = aNorm(Y, V); 

        normVals(i) = aN * (1 + err(2));
        LowBound(i) = normVals(i);

        if prop
            NormErr(i) = dA(1) * normVals(i);
        else
            NormErr(i) = dA(i);
        end

        currEigs = abs(eigs(Y,1,'largestabs'));
        if HighRadius == 0
            HighRadius = currEigs;
        else
            HighRadius = min(HighRadius, currEigs);
        end

        % Accumulate blocks
        X = [X; Y];
    end

    dA        = NormErr;
    NormA     = normVals;
    LowRadius = min(normVals + NormErr);
    mm        = N;   
    iterCount = 1;    
    JJ        = N;    
    MaxJJ     = N;    
    BestPower = 1;
    ell_slp   = 1;

    while (mm < M) && (LowRadius < HighRadius - delta(2))

        iterCount = iterCount + 1;

        OldHighRadius = HighRadius;
        OldLowRadius  = LowRadius;
        LowRadius     = HighRadius;

        NewJJ = 0;
        
        newCapacity        = JJ*N;
        NewNorm            = zeros(newCapacity,1);
        NewNormErr         = zeros(newCapacity,1);
        NewLowBound        = zeros(newCapacity,1);
        NewX               = zeros(n*newCapacity, n);

        newIndex = 0; 

        for k = 1:JJ
            % Extract block #k from X
            Xk = X((k-1)*n+1:k*n, :);

            for i = 1:N
                newIndex = newIndex + 1;

                XX = Xk * A(:,(i-1)*n + 1 : i*n);
                
                NewNormErr(newIndex) = ...
                    ( normVals(k)*(err(1)*NormA(i) + dA(i)) + ...
                      NormErr(k) * (NormA(i) + dA(i)) );

                % antinorm
                [aN,~] = aNorm(XX, V);
                NewNorm(newIndex) = aN * (1 + err(2));

                % new lower bound candidate
                cand = (NewNorm(newIndex) + NewNormErr(newIndex))^(1/iterCount);
                NewLowBound(newIndex) = max(LowBound(k), cand);

                % update HighRadius
                newEig = abs(eigs(XX,1,'largestabs'));
                HighRadius = min(HighRadius, newEig^(1/iterCount));

                if NewLowBound(newIndex) < HighRadius - delta(1)
                    LowRadius = min(LowRadius, NewLowBound(newIndex));
                    NewX((NewJJ)*n + 1 : (NewJJ+1)*n, :) = XX;
                    NewNorm(NewJJ+1)    = NewNorm(newIndex);
                    NewNormErr(NewJJ+1) = NewNormErr(newIndex);
                    NewLowBound(NewJJ+1)= NewLowBound(newIndex);
                    NewJJ = NewJJ + 1;
                end

            end
        end

        X        = NewX(1 : NewJJ*n, :);
        NormA    = NewNorm(1:NewJJ);
        normVals = NormA;
        NormErr  = NewNormErr(1:NewJJ);
        LowBound = NewLowBound(1:NewJJ);

        LowRadius = max(OldLowRadius, min(LowRadius, HighRadius - delta(1)));

        mm = mm + JJ*N;
        JJ = NewJJ;
        MaxJJ = max(JJ, MaxJJ);

        if (HighRadius - LowRadius) < (OldHighRadius - OldLowRadius)
            BestPower = iterCount;
        end

        if HighRadius < OldHighRadius
            ell_slp = iterCount;
        end

        if Display
            disp([BestPower, iterCount, mm, JJ]);
            disp([LowRadius, HighRadius]);
        end

    end

    lsr = [LowRadius, HighRadius];

    m  = [BestPower, iterCount, mm, MaxJJ, ell_slp];


end
