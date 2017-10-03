function [h] = calc_h(R_gas,P,T,Pc,Tc,w,tk,x,ref,mu1,ci2)

N = numel(w);
h = zeros(N-1,1);
mu = zeros(N,1);
lnphi = zeros(N,1);
phi = zeros(N,1);

[lnphi,phi] = Matlab_fugacity(P,T,Pc,Tc,w,tk,x);
[mu] = calc_mu_from_lnphi(R_gas,P,T,x,lnphi);

for j=1:N-1
    [i] = calc_index(ref,j);
    h(j) = ci2(i)*(mu(ref,1) - mu1(ref,1)) - ci2(ref)*(mu(i,1) - mu1(i,1));
end