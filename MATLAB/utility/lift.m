function idx_b = lift(idx_ab, idx_a)
%lift computes indices idx_b such that idx_ab == idx_a(idx_b); i.e., it
%lifts idx_ab to the image of idx_a

[lgc_ab, perm_ab] = unfind(idx_ab, max(idx_a));
[lgc_a, ~, invperm_a] = unfind(idx_a);
assert(all(lgc_a(lgc_ab)), ...
    'Lift not possible if some elements of idx_ab are absent from idx_a')

lgc_b = lgc_ab(lgc_a);
idx_b = find(lgc_b);
idx_b = invperm_a(idx_b(perm_ab));
end