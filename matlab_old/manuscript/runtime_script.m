% This script computes runtime for various matrix operations involving
% LDGMs vs. correlation matrices

clear;clc

% Parallelized over LD blocks
jobnumber = str2num(getenv('SGE_TASK_ID'));
if isempty(jobnumber);jobnumber=1;end


addpath(genpath('..'));

% where to find precision matrices: pmdir/*.EUR.edgelist 
pop = 'EUR';
pm_dir = '/broad/oconnor/trees/nygc/precisionMatrices/';

% where to save results
save_path = '/broad/oconnor/luke/results/070122/matrix_operations_runtime';

% proportion of missing SNPs
missingness = 0.1; 

% replicates
reps = 5;

% these values don't affect runtime
nn = 1e4;
h2 = 0.1;

precision_fnames = dir([pm_dir,'*.',pop,'.edgelist']);

% Load precision matrix data
P = readedgelist([pm_dir,precision_fnames(jobnumber).name]);

% Get rid of missing rows/columns
P = P(any(P), any(P));
noSNPs = length(P);

% Correlation matrix
R = inv(P);

time = zeros(10,2);
for rep = 1:reps

    % Missing data
    notMissing = rand(noSNPs,1) < 1 - missingness;
    
    % Simulate sumstats
    if 1
        tic;
        [sumstats, beta] = simulate_sumstats_populations(...
            nn, {0.5 * ones(length(precisionEstimate),1)},...
            'precisionMatrices',{precisionEstimate},...
            'SNPs',{(1:noSNPs)'},...
            'heritability',h2);
        time(3,2) = time(3,2) + toc;
        
        tic;
        [sumstats0, beta0] = simulate_sumstats_populations(...
            nn, {0.5 * ones(length(precisionEstimate),1)},...
            'correlationMatrices',{R},...
            'SNPs',{(1:noSNPs)'},...
            'heritability',h2);
        time(3,1) = time(3,1) + toc;
                
        sumstats = sumstats{1}(notMissing);
        beta = beta{1};
    end
    
    % Multiply by R
    if 1
        tic;
        x0 = R(notMissing, notMissing) * sumstats;
        time(1,1) = time(1,1) + toc;
        
        % Divide by P
        tic;
        x = precisionDivide(precisionEstimate, sumstats, notMissing);
        time(1,2) = time(1,2) +toc;        
    end
    
    % Divide by R
    if 1
        tic;
        xx = R(notMissing, notMissing) \ sumstats;
        time(2,1) = time(2,1) + toc;
        
        % Multiply by P
        tic;
        x = precisionMultiply(precisionEstimate, sumstats, notMissing);
        time(2,2) = time(2,2) + toc;
    end
    
    % Likelihood given beta using R
    if 1
        tic;
        mu = R(notMissing,:) * beta;
        ll0 = mvnpdf(sumstats, ...
            mu, 1/nn * R(notMissing,notMissing));
          % This approach is slower
%         mu = R(notMissing,:) * beta;
%         x = sumstats0 - mu;
%         ll0 = -1/2 * (log(det(R(notMissing,notMissing))) - log(nn)) - ...
%             1/2 * x' * (R(notMissing,notMissing) \ x);
        time(4,1) = time(4,1) + toc;
        
        % Likelihood given beta using P
        tic;
        mu = precisionEstimate \ beta; mu = mu(notMissing);
        ll = GWASlikelihood(sumstats - mu, zeros(sum(notMissing),1), ...
            precisionEstimate, nn, notMissing);
        time(4,2) = time(4,2) + toc;
    end
    
    % Likelihood given Sigma using R
    if 1
        Sigma = mean(beta.^2) * ones(sum(notMissing),1);
        tic;
        ll0 = mvnpdf(sumstats, zeros(size(sumstats)), ...
            1/nn * R(notMissing,notMissing) + ...
            R(notMissing,notMissing) * diag(Sigma) * R(notMissing,notMissing));
        time(5,1) = time(5,1) + toc;
        
        % Likelihood given Sigma using P
        tic;
        ll = GWASlikelihood(sumstats, Sigma, precisionEstimate, nn, notMissing);
        time(5,2) = time(5,2) + toc;
    end
    
    % Compute the determinant
    if 1
        tic;
        A = chol(R(notMissing,notMissing));
        logDetR0 = 2*sum(log(diag(A)));
        % Slower
%         logDetR0 = log(det(R(notMissing,notMissing)));
        time(6,1) = time(6,1) + toc;
        
        tic;
        A = chol(precisionEstimate);
        B = chol(precisionEstimate(~notMissing,~notMissing));
        logDetR = -2 * (sum(log(diag(A))) - sum(log(diag(B)))) ;
        time(6,2) = time(6,2) + toc;
    end
    
    % BLUP
    if 1
        tau = 1/mean(beta.^2);
        tic;
        mu0 = BLUPslow(sumstats, R, nn, tau, notMissing);
        time(7,1) = time(7,1) + toc;
        
        tic;
        mu = BLUP(sumstats, precisionEstimate, nn, tau, notMissing);
        time(7,2) = time(7,2) + toc;
        
    end
    
    tic;mu = BLUPx({sumstats}, {precisionEstimate}, nn, 1/tau,...
        {1:length(precisionEstimate)}, {notMissing}, {ones(sum(notMissing),1)});
    
    % BLUP variance
    if 1
        tic;
        varBeta0 = diag(inv(nn*R(notMissing,notMissing) + tau.*speye(sum(notMissing))));
        time(8,1) = time(8,1) + toc;
        
        tic;
        denominator = tau * precisionEstimate + nn * speye(noSNPs);
        varBeta = sum(precisionEstimate.*sparseinv(denominator),1)';
        time(8,2) = time(8,2) + toc;
    end
    
    % Sample from posterior
    if 1
        Tau = diag(tau * sparse(notMissing));
        
        tic;
        Tinv = inv(Tau(notMissing,notMissing) + nn*R(notMissing,notMissing));
        x0 = mvnrnd(nn*Tinv * sumstats, Tinv);
        time(9,1) = time(9,1) + toc;
        
        tic;
        x = samplePosteriorMVN(precisionEstimate, chol(precisionEstimate),...
            Tau, nn, sumstats, notMissing);
        time(9,2) = time(9,2) + toc;
    end
    
    % Likelihood gradient
    if 1
        tic;
        Sigma = mean(beta.^2) * ones(sum(notMissing),1);
        RSinv = inv(R(notMissing,notMissing).*Sigma' + 1/nn * eye(sum(notMissing)));
        RSinv_sumstats = RSinv * sumstats;
        grad0 = -1/2 * sum(RSinv .* R(notMissing,notMissing))' + 1/2 * sum(RSinv_sumstats .* RSinv_sumstats')';
        time(10,1) = time(10,1) + toc;
        
        tic;
        grad = GWASlikelihoodGradient(...
            sumstats, 1/tau, precisionEstimate, nn, 1, notMissing, 1);
        time(10,2) = time(10,2) + toc;
    end
    disp(time)
end
time = time/reps;
save([save_path, '.', num2str(jobnumber), '.mat'],'time','noSNPs');













