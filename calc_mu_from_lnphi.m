function [mu] = calc_mu_from_lnphi(R_gas,P,T,x,lnphi)

mu = lnphi + log(P/1e05) + log(x);
mu = R_gas*T*mu;