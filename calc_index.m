function [i] = calc_index(ref,j)

if (j < ref)
    i = j;
else
    i = j+1;
end