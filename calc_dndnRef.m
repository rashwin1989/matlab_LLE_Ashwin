function [dndnRef] = calc_dndnRef(R_gas,P,T,Pc,Tc,w,tk,n,ref,ci2)

N = numel(w);
dndnRef = zeros(N-1,1);

[dmu_dn] = calc_dmu_dn(R_gas,P,T,Pc,Tc,w,tk,n);

A = zeros(N-1,N-1);
b = zeros(N-1,1);

for k=1:N-1    
    for l=1:N-1     
        [i] = calc_index(ref,k);
        [j] = calc_index(ref,l);
        A(k,l) = ci2(i)*dmu_dn(ref,j) - ci2(ref)*dmu_dn(i,j);
        b(k) = ci2(ref)*dmu_dn(i,ref) - ci2(i)*dmu_dn(ref,ref);
    end
end

dndnRef = linsolve(A,b);