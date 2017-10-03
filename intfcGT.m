clear all;
close all;

sp=[11 19];

P0=30*1e06;
T=670;
m12=1;

Nstep=11;
ref=1;
ciType=1;

tol = 1e-8;
iters_max = 108;
R_gas = 8.3144621;

d=load('petro.dat+');
Tc=d(sp,4);
Pc=d(sp,5)*1e5;
w =d(sp,6);
MW=d(sp,7);
for i=1:numel(MW)
 invMW(i,1) = 1/MW(i,1);    
end
N=numel(w);
tk=-sp';

[A,B] = calculate_A_B(w);

if (ciType<0)
    d=load('infl_param.dat+');
    ci=d(:,1);
else    
    [ci]=calculate_ci(T,Pc,Tc,w,A,B);
end
ci2 = ci.^0.5;

d=load('y0.dat+');
y01=d(1,:);
y02=d(2,:);
x01=y01'.*invMW/(y01*invMW);
x02=y02'.*invMW/(y02*invMW);

x1=x01;
x2=x02;

%%Phase-split calculation%%
m1=1/(1 + 1/m12);
m2=1-m1; 
c0=m1/(x1'*MW)/(m1/(x1'*MW) + m2/(x2'*MW));
display([num2str(c0) ': mixing ratio (wt) of crude vs. water is ' ...
         num2str(m12) ':1' ]);

[x1 x2 c]=Matlab_mLLE_new(c0,P0,T,Pc,Tc,w,tk,x01,x02);

[v1]=Matlab_volume(P0,T,Pc,Tc,w,tk,x1);
[v2]=Matlab_volume(P0,T,Pc,Tc,w,tk,x2);
mu1 = zeros(N,1); lnphi = zeros(N,1); phi = zeros(N,1);
[lnphi,phi] = Matlab_fugacity(P0,T,Pc,Tc,w,tk,x1);
[mu1] = calc_mu_from_lnphi(R_gas,P0,T,x1,lnphi);
[lnphi,phi] = Matlab_fugacity(P0,T,Pc,Tc,w,tk,x2);
[mu2] = calc_mu_from_lnphi(R_gas,P0,T,x2,lnphi);

%%Gradient Theory interface calculation%%
n=zeros(N,Nstep);
for i=1:N
    n(i,1)=x2(i)/v2;
    n(i,Nstep)=x1(i)/v1;
end

dnref = ( n(ref,Nstep) - n(ref,1) )/( Nstep - 1 );
P = zeros(Nstep,1);
v = zeros(Nstep,1);
x = zeros(N,Nstep);
P = P0;
dndnRef = zeros(N,1);
dn_dnRef = zeros(N-1,1); 
dnk = zeros(N-1,1);
dzdn = zeros(Nstep,1);
dsigdn = zeros(Nstep,1);
muk = zeros(N,1);
h = zeros(N-1,1);
dhdn = zeros(N-1,N-1);
norm_h = 1.0;

dOmega_const = n(:,Nstep)'*mu1 - P0;

for k=1:2

    nk = n(:,k);
    sum_nk = nk'*ones(N,1);
    vk = 1/sum_nk;
    v(k) = vk;
    xk = nk*vk;
    x(:,k) = xk;
    [Pk] = Matlab_pressure(vk,T,Pc,Tc,w,tk,xk);
    P(k) = Pk;
    [dmudn] = calc_dmu_dn(R_gas,Pk,T,Pc,Tc,w,tk,nk);
    [dn_dnRef] = calc_dndnRef(R_gas,Pk,T,Pc,Tc,w,tk,nk,ref,ci2);
    for j=1:N-1
        [i] = calc_index(ref,j);
        dndnRef(i) = dn_dnRef(j);
    end
    dndnRef(ref) = 1.0;
    dzdnk = 0;
    for i=1:N
        for j=1:N
            dzdnk = dzdnk + 0.5*ci2(i)*ci2(j)*dndnRef(i)*dndnRef(j);
        end
    end
    [lnphi,phi] = Matlab_fugacity(Pk,T,Pc,Tc,w,tk,xk);
    [muk] = calc_mu_from_lnphi(R_gas,Pk,T,xk,lnphi);
    dOmega = nk'*muk - Pk - nk'*mu1 + P0;
    dzdn(k) = (dzdnk/dOmega)^0.5;
    dsigdn(k) = 2.0*dOmega*dzdn(k);

    nkOld = nk;
    nk = nkOld + dndnRef*dnref;

    norm_h = 1.0;
    iters = 0;    
    while ((norm_h > tol) & (iters < iters_max))
        iters = iters + 1;
        nkOld = nk;
        sum_nk = nk'*ones(N,1);
        vk = 1/sum_nk;
        xk = nk*vk;
        [Pk] = Matlab_pressure(vk,T,Pc,Tc,w,tk,xk);
        [h] = calc_h(R_gas,Pk,T,Pc,Tc,w,tk,xk,ref,mu1,ci2);
        [dhdn] = calc_dhdn(R_gas,Pk,T,Pc,Tc,w,tk,nk,ref,ci2);
        dnk = linsolve(dhdn,h);
        [nk] = correct_nk(nkOld,-dnk,ref);
        norm_h = norm(h,Inf);
    end
    if (k < Nstep)
        n(:,k+1) = nk;
    end
end

sigma = 0;
z = zeros(Nstep,1);

for k=1:(Nstep-1)
    sigma = sigma + dnref*0.5*(dsigdn(k) + dsigdn(k+1));
    z(k+1) = z(k) + dnref*0.5*(dzdn(k) + dzdn(k+1));
end

m = numel(z);

% fig1 = figure(1);
% hold on;      
% plot(z,n(1,:),'-k','linewidth',4);
% plot(z(1:10:m),n(1,1:10:m),'sk','markersize',12, ...
%      'markerfacecolor','k','markeredgecolor','k');
% plot(z,n(2,:),'-b','linewidth',4);
% plot(z(1:10:m),n(2,1:10:m),'ob','markersize',12, ...
%      'markerfacecolor','b','markeredgecolor','b');
% ha=xlabel('n (mol/m3)');
% hb=ylabel('z (m)');
% set([ha hb],'fontsize',30,'fontWeight','bold');
% set(gca,'fontsize',30,'fontWeight','bold');
% NumTicks = 6;
% L = get(gca,'YLim');
% set(gca,'YTick',linspace(L(1),L(2),NumTicks));
% box on;
% savefig(fig1,'n_profiles')