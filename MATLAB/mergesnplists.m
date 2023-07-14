function [whichIndices, mergedSumstats, whichSNPs, sumstats_SNPs_in_snplists] = mergesnplists(snplists,sumstats,P,varargin)
%mergesnplists merges a cell array of .snplist tables, one for each LD block,
%with a table of summary statistics or other SNP information.
%
% Input arguments:
% snplists: cell array where each cell contains a snplist table.
%
% sumstats: summary statistics table with mandatory column name
% (non-case-sensitive) SNP or RSID, to be used for merging.
%   Optionally, sumstats can also contain columns named A1/A2
%   or anc_alleles/deriv_alleles. If these columns are
%   specified, they will be merged with the anc_alleles/deriv_alleles columns
%   of snplists.
%
% P: cell array of LDGM precision matrices. Any indices whose
%  corresponding rows/columns in P are empty will be ingored.
%  Can be the empty cell array, in which case all indices are kept.
% 
% Optional input arguments:
% Optional inputs specify the names of various columns in the summary
% statistics table. Only the variant ID column is required; the others arefunction bgenreader
% optional, and will be ignored if they do not match anything in the
% sumstats table
% 
%   columnNameContainingVariantId: column of sumstats table containing SNP 
%   IDs (default: 'SNP')
% 
%   columnNameContainingReferenceAllele, columnNameContainingAlternativeAllele:
%   column names of reference and alternative alleles,
%   required for phasing effect sizes. Any rows with matching SNP IDs but
%   mismatching alleles will be discarded. (default: 'A1', 'A2')
% 
%   columnNameContainingZ: Z scores column name. If not specified, function
%   will attempt to compute Z scores from 'beta' and 'SE' columns. Either
%   way, it will create a column named 'Z_deriv_allele' for the Z score of
%   the derived allele (if ref and alt alleles are found) (default: 'Z')
% 
%   columnNameContainingBeta, columnNameContainingStandardError: if Z score column not
%   found, function will compute Z scores from these columns (default: 
%   'Beta', 'SE')
% 
%   columnNameContainingAlternativeAlleleFrequency: column name containing
%   allele frequency of the *alternative allele*. If specified, function
%   will create a column called AF_deriv_allele with the AF of the derived
%   allele (if ref and alt alleles are found) (default: 'EAF')
% 
%   columnNameContainingPGSWeight: column name containing polygenic score
%   weights. Like Z scores, these are phased and assigned to a new column
%   called 'pgs_weight_deriv_allele'. However, they are handled differently
%   for SNPs that correspond to the same row/column of the LDGM (which are
%   expected to be in near-perfect LD). Instead of picking one of the SNPs
%   arbitrarily, the entry of pgs_weight_deriv_allele will be the sum of
%   the phased PGS weights for all SNPs corresponding to that index.
%   (default: none)
% 
%   retainEquivalentSNPs: option to keep SNPs belonging to the same row/col
%   index, instead of picking a representative
% 
% Ouput arguments:
% whichIndices: cell array of indices, a subset of the 'index' column of each
% snplist. Indexes which row/column of LDGMs have a corresponding SNP in
% the sumstats table. By default, indices are unique, unless
% retainEquivalentSNPs is set to true
%
% mergedSumstats: cell array of sumstats tables, each containing a subset
% of the rows of the original sumstats table, that match the indices of
% each LD block. Note: when multiple SNPs correspond to the same row/column
% of the LDGM, the first-listed one will arbitrarily be chosen as a
% representative, except for the column specified by columnNameContainingPGSWeight
%
% If alleles are specified in the sumstats table, the mergedSumstats tables
% will have an extra column appended called 'phase', which indicates
% whether the alleles match (+1) or anti-match (-1). If they mismatch
% (i.e., flipping their labels doesn't make them match), then the SNP will
% be discarded.
%
% Note: it makes this function run faster if you go
% chromosome-by-chromosome; this is not done automatically.

p = inputParser;

% SNP lists for each LD block
addRequired(p, 'snplists', @iscell);

% summary statistics table
addRequired(p, 'sumstats', @istable);

% cell array of precision matrices or boolean indices indicating which
% indices for each LD block should be retained
addRequired(p, 'P', @iscell);

% optional name of file format. Currently only supported option is 'vcf'
addParameter(p, 'tableFormat', '', @ischar);

% column of sumstats table containing SNP IDs
addParameter(p, 'columnNameContainingVariantId', 'SNP', @(x)ischar(x) || isnumeric(x));

% column of sumstats table containing reference allele
addParameter(p, 'columnNameContainingReferenceAllele', 'A1', @(x)ischar(x) || isnumeric(x));

% column of sumstats table containing alternative allele
addParameter(p, 'columnNameContainingAlternativeAllele', 'A2', @(x)ischar(x) || isnumeric(x));

% column of sumstats table containing Z scores
addParameter(p, 'columnNameContainingZ', 'Z', @(x)ischar(x) || isnumeric(x));

% column of sumstats table containing effect sizes to compute Z scores
addParameter(p, 'columnNameContainingBeta', 'Beta', @(x)ischar(x) || isnumeric(x));

% column of sumstats table containing effect size SE to compute Z scores
addParameter(p, 'columnNameContainingStandardError', 'SE', @(x)ischar(x) || isnumeric(x));

% column of sumstats table containing allele frequencies of alt alleles
addParameter(p, 'columnNameContainingAlternativeAlleleFrequency', 'EAF', @(x)ischar(x) || isnumeric(x));

% whether to match SNPs by position or by variant ID. If true, requires
% that summary statistics contain chromosome and position identifiers, and
% that snplists contain chromosome and position as well (with the same
% genome build).
addParameter(p, 'matchByPosition', false, @islogical);

% column of sumstats table containing chromosome of each SNP
addParameter(p, 'columnNameContainingChromosome', 'CHR', @(x)ischar(x) || isnumeric(x));

% column of sumstats table containing position of each SNP
addParameter(p, 'columnNameContainingPosition', 'POS', @(x)ischar(x) || isnumeric(x));

% optional QC step is to compare allele frequency from sumstats file
% with AF from snplist file, and discard SNPs that seem to have discordant
% frequencies. If doing this, specify which population (e.g., 'EUR') to
% check against in the snplist tables. Optionally specify a tolerance for
% abs(AF1 - AF2).
addParameter(p, 'populationForCheckingAlleleFrequency', '', @(x)ischar(x) || isnumeric(x));
addParameter(p, 'alleleFrequencyDifferenceTolerance', 0.2, @(x)ischar(x) || isnumeric(x));

% columns of sumstats table containing PGS weights
addParameter(p, 'columnNameContainingPGSWeight', '', @(x)ischar(x) || isnumeric(x));

% option to retain equivalent SNPs. Equivalent SNPs are assigned to the
% same row/column of the LDGM, and they are typically in near-perfect LD
% except for imputation issues. With this option set to true, whichIndices
% will have duplicate entries, indicating which rows of mergedSumstats
% belong to the indicated SNPs.
addParameter(p, 'retainEquivalentSNPs', false, @islogical);

parse(p, snplists, sumstats, P, varargin{:});

% turns variables named p.Results.x into just x
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end
clear p

assert(iscell(snplists), 'Please specify snplists as a cell array of tables')
assert(all(cellfun(@istable,snplists)), 'Please specify snplists as a cell array of tables')
assert(istable(sumstats), 'Please specify sumstats as a table')

if strcmpi(tableFormat, 'vcf')
    columnNameContainingVariantId = 'ID';
    columnNameContainingReferenceAllele = 'A1';
    columnNameContainingAlternativeAllele = 'A2';
    columnNameContainingBeta = 'BETA';
    columnNameContainingStandardError = 'SE';
end

% Identify columns of sumstats table
sumstatsColumnNames = sumstats.Properties.VariableNames;
snpcolumn = getMatchingColumn(sumstatsColumnNames,columnNameContainingVariantId);
poscolumn = getMatchingColumn(sumstatsColumnNames,columnNameContainingPosition);
chrcolumn = getMatchingColumn(sumstatsColumnNames,columnNameContainingChromosome);
zcolumn = getMatchingColumn(sumstatsColumnNames, columnNameContainingZ);
a1column = getMatchingColumn(sumstatsColumnNames, columnNameContainingReferenceAllele);
a2column = getMatchingColumn(sumstatsColumnNames, columnNameContainingAlternativeAllele);
betacolumn = getMatchingColumn(sumstatsColumnNames, columnNameContainingBeta);
secolumn = getMatchingColumn(sumstatsColumnNames, columnNameContainingStandardError);
afcolumn = getMatchingColumn(sumstatsColumnNames, columnNameContainingAlternativeAlleleFrequency);
pgscolumn = getMatchingColumn(sumstatsColumnNames, columnNameContainingPGSWeight);

if matchByPosition
    assert(sum(poscolumn) == 1, ...
        'Please specify exactly one column in sumstats table containing variant positions')
    assert(sum(chrcolumn) == 1, ...
        'Please specify exactly one column in sumstats table containing chromosome of each variant')
else
    assert(sum(snpcolumn) == 1, ...
        'Please specify exactly one column in sumstats table containing variant identifiers')
end
concatenated_snplists = vertcat(snplists{:});

if matchByPosition
    [~, ldgm_idx, sumstats_idx] = intersect([concatenated_snplists.chr, concatenated_snplists.position],...
        [table2array(sumstats(:,chrcolumn)) table2array(sumstats(:,poscolumn))], 'rows', 'stable');
else
    [~, ldgm_idx, sumstats_idx] = intersect(concatenated_snplists.site_ids,...
        table2cell(sumstats(:,snpcolumn)), 'stable');
end

% Subset sumstats to matching SNPs
sumstats = sumstats(sumstats_idx,:);

% Indices to recover LD blocks from concatenated list
blocksizes = cellfun(@height,snplists);
cumulative_blocksizes = [0;cumsum(blocksizes)];
blocks = arrayfun(@(i,j){i+1:j}, cumulative_blocksizes(1:end-1), cumulative_blocksizes(2:end));

% Which SNPs in the snplists have a corresponding entry in the
% sumstats file
whichSNPs = false(cumulative_blocksizes(end),1);
whichSNPs(ldgm_idx) = true;
whichSNPs = cellfun(@(i){whichSNPs(i)}, blocks);

% Indices to recover merged LD blocks from concatenated list
blocksizes = cellfun(@sum,whichSNPs);
cumulative_blocksizes = [0;cumsum(blocksizes)];
blocks = arrayfun(@(i,j){i+1:j}, cumulative_blocksizes(1:end-1), cumulative_blocksizes(2:end));

% Sumstats tables merged with each snplist
mergedSumstats = cellfun(@(ii){sumstats(ii,:)}, blocks);

% which sumstats SNPs had a matching SNP in each LD block
sumstats_SNPs_in_snplists = cellfun(@(i){sumstats_idx(i)}, blocks);


% Phase alleles between snplists and sumstats
noBlocks = numel(snplists);
if sum(a1column)==1 && sum(a2column) == 1
    for ii = 1:noBlocks

        idx = find(whichSNPs{ii});

        % +1 for matching alleles, -1 for anti-matching, 0 for mismatching
        phase = mergealleles(table2cell(mergedSumstats{ii}(:,a1column)), ...
            table2cell(mergedSumstats{ii}(:,a2column)), ...
            snplists{ii}.anc_alleles(idx),...
            snplists{ii}.deriv_alleles(idx));
        mergedSumstats{ii}.phase = phase;

        % get rid of mismatched alleles
        if mean(phase==0) > 0.5
            warning('In block %d, more than half of putatively matching SNPs had mismatched alleles (perhaps due to strandedness?)',ii)
        end
        mergedSumstats{ii} = mergedSumstats{ii}(phase~=0, :);
        whichSNPs{ii} = idx(phase~=0, :);
    end
else
    for ii = 1:noBlocks
        mergedSumstats{ii}.phase = ones(height(mergedSumstats{ii}),1);
    end
    warning('Did not find allele information in sumstats file. Please ensure that alleles are matched between summary statistics and the LDGM')
end


% Compute derived allele Z scores, derived AFs, derived allele PGS weights
for ii = 1:noBlocks
        phase = mergedSumstats{ii}.phase;
        
        % Z score of derived allele
        if any(zcolumn)
            mergedSumstats{ii}.Z_deriv_allele = phase .* ...
                table2array(mergedSumstats{ii}(:,zcolumn));
        elseif any(betacolumn) && any(secolumn)
            mergedSumstats{ii}.Z_deriv_allele = phase .* ...
                table2array(mergedSumstats{ii}(:,betacolumn)) ./...
                table2array(mergedSumstats{ii}(:,secolumn));
        end

        % AF of derived allele
        if any(afcolumn)
            mergedSumstats{ii}.AF_deriv_allele = (1-phase)/2 + phase .* ...
                table2array(mergedSumstats{ii}(:,afcolumn));
        end

        % PGS weight of derived allele
        if any(pgscolumn)
            mergedSumstats{ii}.pgs_weight_deriv_allele = phase .* ...
                table2array(mergedSumstats{ii}(:,pgscolumn));
        end
end

% Convert from SNPs to LDGM row/col indices (which can have multiple SNPs)
whichIndices = cell(size(snplists));
if retainEquivalentSNPs
    for ii = 1:noBlocks
        whichIndices{ii} = snplists{ii}.index(whichSNPs{ii}) + 1;

        % Get rid of indices whose corresponding columns of P are empty
        if ~isempty(P)
            include = ismember(whichIndices{ii}, find(any(P{ii})));
            whichIndices{ii} = whichIndices{ii}(include);
            mergedSumstats{ii} = mergedSumstats{ii}(include,:);
        end
    end
else
    representatives = cell(size(snplists));
    for ii = 1:noBlocks
        [whichIndices{ii}, representatives{ii}, duplicates] = ...
            unique(snplists{ii}.index(whichSNPs{ii}) + 1);

        % For PGS column, add up the entries when multiple SNPs are
        % assigned to the same index
        if any(pgscolumn)
            mergedSumstats{ii}.pgs_weight_deriv_allele(representatives{ii},:) = ...
                accumarray(duplicates, ...
                mergedSumstats{ii}.pgs_weight_deriv_allele);
        end

        % Get rid of indices whose corresponding columns of P are empty
        if ~isempty(P)
            include_indices = find(any(P{ii}));
            [whichIndices{ii}, include_indices] = intersect(whichIndices{ii},include_indices);
            representatives{ii} = representatives{ii}(include_indices);
        end

        mergedSumstats{ii} = mergedSumstats{ii}(representatives{ii},:);
    end
    
    % Get rid of SNPs whose allele frequency is discordant
    if ~isempty(populationForCheckingAlleleFrequency)
        for block = 1:noBlocks
            AF1 = mergedSumstats{block}.AF_deriv_allele;
            AF2 = table2array(snplists{block}(whichSNPs{block}, populationForCheckingAlleleFrequency));
            AF2 = AF2(representatives{block});
            AF_difference = abs(AF1 - AF2);
            incl = AF_difference < alleleFrequencyDifferenceTolerance;
            if mean(incl) < .9
                warning('Large number of SNPs (>0.1) with discordant allele frequencies')
            end

            whichIndices{block} = whichIndices{block}(incl);
            mergedSumstats{block} = mergedSumstats{block}(incl,:);
            whichSNPs{block} = whichSNPs{block}(incl);

        end
    end
end

end

function columnMatch = getMatchingColumn(columnNames, columnSpecifier)
    if isempty(columnSpecifier)
        columnMatch = false;
        return;
    end
    if isnumeric(columnSpecifier)
        assert(all(round(columnSpecifier) == columnSpecifier), ...
            'Column specifiers should be strings or indices')
        assert(all(columnSpecifier >= 1), ...
            'Column specifiers should be strings or indices')
        assert(all(columnSpecifier <= length(columnNames)), ...
            'Column specifier indices should be at most the number of columns in the summary statistics table')
        columnMatch = unfind(columnSpecifier, length(columnNames));
        return;
    end
    if ischar(columnSpecifier)
        columnMatch = strcmpi(columnSpecifier, columnNames);
        return;
    end
    error('Column specifiers should be strings or indices')
    
end

