function [SNPs,chisq,beta,se,freq,info] = import_sumstats_boltfile(filepath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
file=fopen(filepath);
fgets(file);

% you might have to change this line (eg because there is no INFO column)
data=textscan(file,'%s %*d %*d %*d %s %s %f %f %*f %*f %f %f %*f %*f %f %*f');

if ~feof(file)
    error('Failed to load summary stats')
end
fclose(file);
freq=data{4};
info=data{5};
beta=data{6};
se=data{7};
chisq=data{8};

% Truncates "rs" and returns 0 for any SNPs that lack RSIDs
SNPs=cellfun(@(s)scan_rsid(s),data{1});

% aligns beta to common reference allele using alphabetical order of
% allele1 and allele2. This avoids issues between datasets where sometimes
% the alternative/effect allele will be set to the minor allele, which can
% differ across datsets when the AF is very close to 50%
phase=(-1).^(cellfun(@(x,y)mod(sum(x(1)<y(1)),2)==0,data{2},data{3}));

beta=beta.*phase;

    function number=scan_rsid(string)
        number=sscanf(string,'rs%d');
        if isempty(number)
            number=0;
        end
        
    end
end

