function [betaPosteriorMean] =...
    gibbsGWAS(P, alphaHat, nn, tauSupport, tauPrior, reps, blocks)
%betaExpectation: expected value of beta given the data
% alphahat: GWAS sumstats; omega: LD precision matrix; nn: GWAS sample
% size; tauSupport: possible values for tau; tauPrior: prior probability of
% each tau value; reps: number of steps to take; betasq: don't use it

% betaHat = P * alphaHat;
mm = length(P);
cholP = chol(P)';
if nargin < 8
    blocks = {1:mm};
end
noBlocks = length(blocks);

if sum(tauPrior)~=1
    error('Prior should sum to 1')
end

mm = length(alphaHat);

% Initialize
tauBLUP = 1/sum(tauPrior./tauSupport);
beta = BLUP(alphaHat, P, nn, tauBLUP);

betaPosteriorMean = zeros(mm,1);
betaMean = betaPosteriorMean;
for rep = 1 : reps
    likelihood = sqrt(tauSupport'/(2*pi)) .* ...
        exp(-1/2 * tauSupport' .* beta.^2);
    
    tauPosterior = tauPrior' .* likelihood;
    tauPosterior = tauPosterior ./ sum(tauPosterior, 2);
    
    [~, whichTau] = max(rand(mm,1) < ...
        cumsum(tauPosterior, 2), [], 2);
    tau = tauSupport(whichTau);
    
    for b = 1:noBlocks
        [beta(blocks{b}), betaMean(blocks{b})] = ...
            samplePosteriorMVN(P(blocks{b},blocks{b}),...
            cholP(blocks{b},blocks{b}), speye(length(blocks{b})).*tau(blocks{b}),...
            nn, alphaHat(blocks{b}));
    end
    %     betaPosteriorVariance = diag(1./( tau + nn));
    %     betaMean = nn * betaPosteriorVariance * alphaHat;
    %     beta = sqrt(betaPosteriorVariance) * randn(mm,1) + betaMean;
    
    betaPosteriorMean = betaPosteriorMean + betaMean;
    
    
end

betaPosteriorMean = betaPosteriorMean / rep;
end

