function [lsr, m, err] = lsr_comp(A, dA, delta, M, V, err, Display)
%
%    Inputs (with default values):
%       A       : (n x n*N) or (n*N x n) matrix family elements (mandatory).
%       dA      : estimate of the norm of the error in A's matrices (default []).
%       delta   : two-element vector controlling tolerance [delta_low, delta_high] (default 1e-2).
%       M       : max number of expansions (default 250).
%       V       : matrix for aNorm computations (default identity).
%       err     : relative error array for aNorm computations (default [2^(-50)*n, 2^(-49)*n]).
%       rho     : initial guess for the largest stable radius? (default 0).
%       Display : whether to display iteration info (default 0).
%
%    Outputs:
%       lsr : [LowRadius, HighRadius] final bounds.
%       m   : [BestPower, m_current, mm, MaxJJ, ell_slp] iteration info.

    if nargin < 2 || isempty(dA),   dA = [];           end
    if nargin < 3 || isempty(delta), delta = 1e-2;     end
    if nargin < 4 || isempty(M),     M = 250;          end
    if nargin < 5 || isempty(V),     V = eye(size(A,1));  end
    if nargin < 6 || isempty(err),   err = [];         end
    if nargin < 7 || isempty(Display), Display = 0;    end

    delta = delta(:);
    if isscalar(delta)
       delta(2) = delta(1);
    end

    [ma, na] = size(A);
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

    if numel(err) < 2
        err = [2^(-50)*n, 2^(-49)*n]; 
    end

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
