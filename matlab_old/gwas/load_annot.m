function [annot,header,SNPs] = load_annot(filename,thin_annot)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if nargout>2
    if thin_annot
        error('Cannot output SNP names from thin annot file')
    end
end
file=fopen(filename);
header=strsplit(fgets(file));
if thin_annot
    header=header(1:end-1);
    annot=textscan(file,[repmat('%f ',1,length(header)), '\n']);
    annot=[annot{1:end}];
else
    header=header(5:end);
    annot=textscan(file,['%*d %*d %*c%*c%d %*f',repmat('%f ',1,length(header)), '\n']);
    SNPs=annot{1};
annot=[annot{2:end-2}];
end

fclose(file);
end

