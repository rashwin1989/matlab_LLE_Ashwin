function [n] = calc_nk_next(nOld,dndnRef,dnref,ref)

n = nOld;

for j=1:N-1
    [i] = calc_index(ref,j);
    n(i) = n(i) + dndnRef(j)*dnref;
end