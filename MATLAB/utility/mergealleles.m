function phase = mergealleles(A1,A2,B1,B2)
%mergealleles compares alleles A1, A2 with B1, B2 and returns either:
% 1, wherever A1==B1 and A2==B2
% -1, wherever A1==B2 and A2==B1
% 0, otherwise
% All four input arguments should have the same size; comparisons are
% non-case sensitive

assert(all(size(A1)==size(A2)))
assert(all(size(A1)==size(B1)))
assert(all(size(A1)==size(B2)))

% SNPs whose A1/A2 alleles both match
matching_alleles = strcmpi(A1, B1) & ...
    strcmpi(A2, B2);

% SNPs whose A1/A2 alleles crossmatch
flipped_alleles = strcmpi(A1, B2) & ...
    strcmpi(A2, B1);

% SNPs with mismatch
wrong_alleles = ~(matching_alleles | flipped_alleles);

assert(all(matching_alleles + flipped_alleles + wrong_alleles == 1))

phase = zeros(size(wrong_alleles));
phase(matching_alleles) = 1;
phase(flipped_alleles) = -1;

end

