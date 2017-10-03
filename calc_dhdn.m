function [dhdn] = calc_dhdn(R_gas,P,T,Pc,Tc,w,tk,n,ref,ci2)

N = numel(w);
dhdn = zeros(N-1,N-1);
dmudn = zeros(N,N);
[dmudn] = calc_dmu_dn(R_gas,P,T,Pc,Tc,w,tk,n);

for i=1:N-1
    for j=1:N-1
        [k] = calc_index(ref,i);
        [l] = calc_index(ref,j);
        dhdn(i,j) = ci2(k)*dmudn(ref,l) - ci2(ref)*dmudn(k,l);
    end
end