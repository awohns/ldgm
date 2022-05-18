function [l2] = getLDSC(P,A)
%getLDSC compute stratified LD scores from precision matrices and
%annotation matrices
% P: cell array of precision matrices
% A: cell array of annotation matrices; each entry of A should have number
% of rows == number of rows+columns in each entry of P
% l2: LD scores

noAnnot = size(A{1},2);
l2 = cell(size(P));
for b=1:length(P)
    R = inv(P{b});
    for k = 1:noAnnot
        l2{b}(:,k) = sum(R.^2.*A{b}(:,k));
    end

end
end

