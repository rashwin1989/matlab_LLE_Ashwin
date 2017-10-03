function [dmu_dn] = calc_dmu_dn(R_gas,P0,T,Pc,Tc,w,tk,n0)

N = numel(w);
mu1 = zeros(N,1);
mu2 = zeros(N,1);
lnphi = zeros(N,1);
phi = zeros(N,1);
dmu_dn = zeros(N,N);

P=P0; n=n0;
sum_n = n'*ones(N,1);
v = 1/sum_n;
x = n*v;
[lnphi,phi] = Matlab_fugacity(P,T,Pc,Tc,w,tk,x);
[mu1] = calc_mu_from_lnphi(R_gas,P,T,x,lnphi);

for j=1:N             
    n = n0;
    n(j) = n0(j) + 1e-6;
    sum_n = n'*ones(N,1);    
    v = 1/sum_n;
    x = n*v;
    [P] = Matlab_pressure(v,T,Pc,Tc,w,tk,x);
    [lnphi,phi] = Matlab_fugacity(P,T,Pc,Tc,w,tk,x);
    [mu2] = calc_mu_from_lnphi(R_gas,P,T,x,lnphi);

    dmu_dn(:,j) = (mu2 - mu1)/(1e-6);
end