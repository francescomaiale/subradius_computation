function tests = test_adaptive_subradius
% Numerical enclosure regressions; finite budgets need not reach delta.
% Run with: results = runtests('tests/test_adaptive_subradius.m');
%           assertSuccess(results)
% SUBRADIUS_TEST_SOURCE can point to an unmodified checkout for comparison.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
originalPath = path;
originalFolder = pwd;
originalRng = rng;
originalWarnings = warning;
testCase.addTeardown(@() path(originalPath));
testCase.addTeardown(@() cd(originalFolder));
testCase.addTeardown(@() rng(originalRng));
testCase.addTeardown(@() warning(originalWarnings));

sourceRoot = getenv('SUBRADIUS_TEST_SOURCE');
if isempty(sourceRoot)
    sourceRoot = fileparts(fileparts(mfilename('fullpath')));
end
if ~isfolder(sourceRoot)
    error('SubradiusTests:MissingSource', ...
        'Cannot find the requested source folder: %s', sourceRoot);
end
addpath(sourceRoot, '-begin');
cd(sourceRoot);
clear adaptive_subradius_comp adaptive_eigenvectors_subradius_comp aNorm real_antinorm
required = {'adaptive_subradius_comp', ...
    'adaptive_eigenvectors_subradius_comp', 'aNorm', 'real_antinorm'};
for index = 1:numel(required)
    if isempty(which(required{index}))
        error('SubradiusTests:MissingSource', ...
            'Required repository function is missing: %s', required{index});
    end
end
if isempty(which('linprog')) || ~license('test', 'Optimization_Toolbox')
    error('SubradiusTests:MissingOptimizationToolbox', ...
        'These numerical regressions require MATLAB and a licensed Optimization Toolbox (linprog).');
end
% A failed LP must fail the regression, rather than silently supply Inf/zero.
warning('error', 'aNorm:SolverError');
warning('error', 'linprog:NoConvergence');
warning('error', 'linprog:Exitflag');
testCase.TestData.methods = {'A', 'E'};
end

function testPaperReferenceWithAndWithoutPruning(testCase)
blocks = paperFamily();
for method = testCase.TestData.methods
    for prune = [false, true]
        result = runMethod(method{1}, blocks, 100, prune);
        checkKnownValue(testCase, result, 1, ...
            sprintf('paper %s, pruning %d', method{1}, prune));
    end
end
end

function testPaperMatrixOrder(testCase)
% Reordering the menu preserves all product words, but changes queue order.
blocks = fliplr(paperFamily());
for method = testCase.TestData.methods
    for prune = [false, true]
        result = runMethod(method{1}, blocks, 100, prune);
        checkKnownValue(testCase, result, 1, ...
            sprintf('swapped paper %s, pruning %d', method{1}, prune));
    end
end
end

function testPaperScaling(testCase)
% Scaling every generator scales the true subradius by the same factor.
base = paperFamily();
for factor = [0.5, 2]
    blocks = cellfun(@(matrix) factor * matrix, base, 'UniformOutput', false);
    for method = testCase.TestData.methods
        result = runMethod(method{1}, blocks, 100, true);
        checkKnownValue(testCase, result, factor, ...
            sprintf('scaled paper %s, factor %g', method{1}, factor));
    end
end
end

function testStrictlyPositiveKnownFamilies(testCase)
% Ai=ri*Si with positive row-stochastic Si: every product has infinity
% norm and spectral radius equal to the product of its ri, so LSR=min(ri).
% Unequal column sums make the initial 1-antinorm bound non-exact.
families = {
    {0.9 * [0.9, 0.1; 0.7, 0.3], 1.2 * [0.2, 0.8; 0.4, 0.6]}, ...
    {1.1 * [0.8, 0.2; 0.3, 0.7], 0.8 * [0.1, 0.9; 0.6, 0.4], ...
        1.3 * [0.45, 0.55; 0.85, 0.15]}
};
expected = [0.9, 0.8];
for familyIndex = 1:numel(families)
    for method = testCase.TestData.methods
        for prune = [false, true]
            result = runMethod(method{1}, families{familyIndex}, 80, prune);
            checkKnownValue(testCase, result, expected(familyIndex), ...
                sprintf('positive family %d, %s, pruning %d', ...
                    familyIndex, method{1}, prune));
        end
    end
end
end

function testIndependentFiniteProductEnvelopes(testCase)
% Enumerate ALL short products independently, without adaptive vertices,
% retained-product queues, or the repository's antinorm/eigenvalue helpers.
% For nonnegative P, min(sum(P,1)) is its induced 1-antinorm. At every
% degree, its minimum root bounds the LSR below, while spectral roots
% bound it above. A returned enclosure cannot cross either outer bound.
families = {
    {[0.7, 0.2; 0.4, 1.1], [1.0, 0.5; 0.1, 0.4]}, ...
    {[1.2, 0.1; 0.3, 0.7], [0.4, 0.9; 0.8, 0.2], ...
        [0.65, 0.35; 0.25, 0.95]}
};
for familyIndex = 1:numel(families)
    envelope = finiteProductEnvelope(families{familyIndex}, 6);
    for method = testCase.TestData.methods
        result = runMethod(method{1}, families{familyIndex}, 100, true);
        label = sprintf('finite-product family %d, %s', familyIndex, method{1});
        checkOutput(testCase, result, label);
        tolerance = numericalTolerance(envelope(2));
        verifyLessThanOrEqual(testCase, result.bounds(1), envelope(2) + tolerance, ...
            [label ': lower estimate exceeds an independently enumerated spectral upper bound']);
        verifyGreaterThanOrEqual(testCase, result.bounds(2), envelope(1) - tolerance, ...
            [label ': upper estimate is below the independently enumerated column-sum lower bound']);
    end
end
end

function testCoarseToleranceRetainsPreviousLowerBound(testCase)
% The first level enters the search, but all children are discarded at the
% deliberately coarse tolerance. The previous 1-antinorm lower certificate
% must survive; merely clamping a reversed interval to its upper end would
% claim unsupported extra accuracy and would fail this independent value.
blocks = {[2, 0.1; 0.1, 0.5], [0.5, 0.1; 0.1, 2]};
initialLower = min(cellfun(@(matrix) min(sum(matrix, 1)), blocks));
productUpper = max(abs(eig(blocks{1} * blocks{2})))^(1/2);
for method = testCase.TestData.methods
    for prune = [false, true]
        result = runMethod(method{1}, blocks, 100, prune, 1.3);
        label = sprintf('coarse tolerance %s, pruning %d', method{1}, prune);
        checkOutput(testCase, result, label);
        verifyEqual(testCase, result.bounds(1), initialLower, ...
            'AbsTol', numericalTolerance(initialLower), ...
            [label ': retain the independent previous lower certificate']);
        verifyEqual(testCase, result.bounds(2), productUpper, ...
            'AbsTol', numericalTolerance(productUpper), ...
            [label ': retain the independent two-product spectral upper estimate']);
    end
end
end

function blocks = paperFamily()
% Transposed Section 5.1 example, as used in illustrative_example.m.
A1 = [7, 2; 0, 3];
A2 = [2, 0; 4, 8];
product = A1 * A2 * (A1^2 * A2)^2;
reference = max(abs(eig(product)))^(1/8);
blocks = {A1 / reference, A2 / reference};
end

function result = runMethod(method, blocks, budget, prune, delta)
% Reset for every run: eigs can otherwise consume a shared random stream.
rng(20261018, 'twister');
dimension = size(blocks{1}, 1);
packed = [blocks{:}];
if nargin < 5, delta = 1e-6; end
if strcmp(method, 'A')
    [bounds, metrics, ~, vertices, witness] = adaptive_subradius_comp( ...
        packed, 0, delta, budget, eye(dimension), [0; 0], false, prune);
else
    [bounds, metrics, ~, vertices, witness] = adaptive_eigenvectors_subradius_comp( ...
        packed, 0, delta, budget, eye(dimension), [0; 0], false, 1.1, prune);
end
result = struct('bounds', bounds, 'metrics', metrics, ...
    'vertices', vertices, 'witness', witness, 'dimension', dimension);
end

function checkKnownValue(testCase, result, expected, label)
checkOutput(testCase, result, label);
tolerance = numericalTolerance(expected);
verifyLessThanOrEqual(testCase, result.bounds(1), expected + tolerance, ...
    [label ': lower estimate excludes the known subradius']);
verifyGreaterThanOrEqual(testCase, result.bounds(2), expected - tolerance, ...
    [label ': upper estimate excludes the known subradius']);
end

function checkOutput(testCase, result, label)
assertSize(testCase, result.bounds, [1, 2], label);
assertTrue(testCase, all(isfinite(result.bounds)), label);
verifyGreaterThanOrEqual(testCase, result.bounds(1), 0, label);
tolerance = numericalTolerance(result.bounds(2));
verifyLessThanOrEqual(testCase, result.bounds(1), result.bounds(2) + tolerance, ...
    [label ': estimates are reversed']);
assertNumElements(testCase, result.metrics, 5, label);
verifyTrue(testCase, all(isfinite(result.metrics)) && all(result.metrics >= 1), label);
assertSize(testCase, result.witness, [result.dimension, result.dimension], label);
verifyTrue(testCase, all(isfinite(result.witness), 'all'), label);
verifyTrue(testCase, all(isfinite(result.vertices), 'all') && ...
    all(result.vertices >= -tolerance, 'all'), label);
degree = result.metrics(5);
assertTrue(testCase, isfinite(degree) && degree >= 1 && degree == fix(degree), label);
witnessBound = max(abs(eig(result.witness)))^(1/degree);
verifyEqual(testCase, witnessBound, result.bounds(2), 'AbsTol', tolerance, ...
    [label ': X_opt and its reported degree do not reproduce the upper estimate']);
end

function bounds = finiteProductEnvelope(blocks, maxDegree)
products = {eye(size(blocks{1}, 1))};
lower = 0;
upper = Inf;
for degree = 1:maxDegree
    next = cell(1, numel(products) * numel(blocks));
    columnMinimum = Inf;
    spectralMinimum = Inf;
    index = 0;
    for parent = 1:numel(products)
        for generator = 1:numel(blocks)
            index = index + 1;
            product = products{parent} * blocks{generator};
            next{index} = product;
            columnMinimum = min(columnMinimum, min(sum(product, 1))^(1/degree));
            spectralMinimum = min(spectralMinimum, max(abs(eig(product)))^(1/degree));
        end
    end
    lower = max(lower, columnMinimum);
    upper = min(upper, spectralMinimum);
    products = next;
end
bounds = [lower, upper];
end

function tolerance = numericalTolerance(value)
% LP reciprocals use a 1e-10 safeguard; leave room for accumulated roundoff,
% while remaining well below the observed approximately 2.5e-5 regression.
tolerance = 2e-8 * max(1, abs(value));
end
