function [l2, RDiag] = getLDSC(P,a,blocks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
noAnnot = size(a,2);
l2 = zeros(size(a));
RDiag = zeros(size(a,1),1);
for b=1:length(blocks)
    R = inv(P(blocks{b},blocks{b}));
    for k = 1:noAnnot
%         A = speye(length(blocks{b})).*a(blocks{b},k);
%         L = P(blocks{b},blocks{b}) \ A / P(blocks{b},blocks{b});
%         l2(blocks{b},k) = diag(L);
        l2(blocks{b},k) = sum(R.^2.*a(blocks{b},k));
    end
    RDiag(blocks{b}) = diag(R);

end
end

