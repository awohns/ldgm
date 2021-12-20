function [alphaHat,beta,alpha,sigmasq] = ...
    simulate_sumstats_blocks(P, nn, sigmasqSupport, varargin)
    

p=inputParser;
mm = length(P);
addRequired(p, 'P', @iscell)
addRequired(p, 'nn', @isscalar)
addRequired(p, 'sigmasqSupport', @isnumeric)
addOptional(p, 'sigmasqPrior', 1, @(x)all(sum(x,2)==1))
addOptional(p, 'annot', cellfun(@(x){ones(length(x),1)},P), @iscell)
addOptional(p, 'linkFn', @(x)max(x,0), @(f)isa(f,'function_handle'))

parse(p,P, nn, sigmasqSupport,varargin{:});

% turns p.Results.x into just x
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end

noBlocks = length(P);
alphaHat = cell(noBlocks,1);
beta = alphaHat; alpha = alphaHat; sigmasq = alphaHat;

for b = 1:noBlocks
    [alphaHat{b}, beta{b}, alpha{b}, sigmasq{b}] = ...
        simulate_sumstats_precision(P{b}, nn, sigmasqSupport, ...
        'sigmasqPrior',sigmasqPrior, 'annot',annot{b}, 'linkFn',linkFn);
end


end