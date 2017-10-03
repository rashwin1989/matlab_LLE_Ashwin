function [nk] = correct_nk(nkOld,dnk,ref);

N = numel(nkOld);
nk = nkOld;

for j=1:N-1
    [i] = calc_index(ref,j);
    nk(i) = nk(i) + dnk(j);
end