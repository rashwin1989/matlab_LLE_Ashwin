clear all;
close all;

sp=[8 19];

%specific gravity of Shanghai oil: 0.9724
%20g oil in batch reactor, 20ml volume of oil, 80ml 8.5MPa N2
%charged = 7.6g N2; oil:N2 = 2.6

d=load('petro.dat+');
Tb=d(sp,3);
Tc=d(sp,4);
Pc=d(sp,5)*1e5;
w =d(sp,6);
MW=d(sp,7);
Vc=d(sp,8);
SG=d(sp,9);
H8=d(sp,10);
for i=1:numel(MW)
 invMW(i,1) = 1/MW(i,1);    
end

n=numel(w);

d=load('groups.dat+');
dm=d(sp,23);
kappa=zeros(n,1);
kappa(n,1) = 0.076;

tk=-sp';

%d=load('Diaz.dat');
ncomp = 3;
d=load('y0.dat+');
for i=1:ncomp
    y0(i,:)=d(i,:);
    x0(:,i)=y0(i,:)'.*invMW/(y0(i,:)*invMW);
end

Dij_tmp=zeros(n,n);

P=27*1e6;

T=[600:1:720];

for j=1:ncomp
    x = x0(:,j);
    for i=1:numel(T)        

        [v]=Matlab_volume(P,T(i),Pc,Tc,w,tk,x);    

        [mu_tmp lambda_tmp]=Matlab_vis_cond(P,T(i),Pc,Tc,Vc,w, ...
                                            MW,kappa,dm,x,v,Tb,SG,H8); 
        [Dij_tmp]=Matlab_Dij_tlsm_wk(P,T(i),Pc,Tc,Vc,w,tk,MW,x);

        for k=1:n
            Di_tmp(k) = 0;
            for l=1:n
                if (k ~= l)
                    Di_tmp(k) = Di_tmp(k) + x(l)/Dij_tmp(k,l);
                end
            end
            Di_tmp(k) = (1 - x(k))/Di_tmp(k);
        end

        rho(i,j) = 1e-3*(x'*MW)/v;
        mu(i,j) = mu_tmp;
        lambda(i,j) = lambda_tmp;
        Di(i,:,j) = Di_tmp;    
    end
end

m=numel(T);

switch n      

  case 2

    fig1 = figure(1);
    hold on;
    plot(T,rho(:,1),'-r','linewidth',4);
    plot(T(1:10:m),rho(1:10:m,1),'sr','markersize',12, ...
         'markerfacecolor','r','markeredgecolor','r');
    plot(T,rho(:,2),'-g','linewidth',4);
    plot(T(1:10:m),rho(1:10:m,2),'dg','markersize',12, ...
         'markerfacecolor','g','markeredgecolor','g');
    plot(T,rho(:,3),'-b','linewidth',4);
    plot(T(1:10:m),rho(1:10:m,3),'ob','markersize',12, ...
         'markerfacecolor','b','markeredgecolor','b');
    ha=xlabel('T (K)');
    hb=ylabel('rho (kg/m3)');
    set([ha hb],'fontsize',30,'fontWeight','bold');
    set(gca,'fontsize',30,'fontWeight','bold');
    NumTicks = 6;
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks));
    box on;
    savefig(fig1,'rho')

    fig2 = figure(2);
    hold on;
    plot(T,Di(:,1,1),'-r','linewidth',4);
    plot(T(1:10:m),Di(1:10:m,1,1),'sr','markersize',12, ...
         'markerfacecolor','r','markeredgecolor','r');    
    plot(T,Di(:,1,2),'-g','linewidth',4);
    plot(T(1:10:m),Di(1:10:m,1,2),'dg','markersize',12, ...
         'markerfacecolor','g','markeredgecolor','g');    
    plot(T,Di(:,1,3),'-b','linewidth',4);
    plot(T(1:10:m),Di(1:10:m,1,3),'ob','markersize',12, ...
         'markerfacecolor','b','markeredgecolor','b');    
    ha=xlabel('T (K)');
    hb=ylabel('D (m2/s)');
    set([ha hb],'fontsize',30,'fontWeight','bold');
    set(gca,'fontsize',30,'fontWeight','bold');
    NumTicks = 6;
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks));
    box on;
    savefig(fig2,'D1')

    fig3 = figure(3);
    hold on;
    plot(T,Di(:,n,1),'-r','linewidth',4);
    plot(T(1:10:m),Di(1:10:m,n,1),'sr','markersize',12, ...
         'markerfacecolor','r','markeredgecolor','r');    
    plot(T,Di(:,n,2),'-g','linewidth',4);
    plot(T(1:10:m),Di(1:10:m,n,2),'dg','markersize',12, ...
         'markerfacecolor','g','markeredgecolor','g');    
        plot(T,Di(:,n,3),'-b','linewidth',4);
    plot(T(1:10:m),Di(1:10:m,n,3),'ob','markersize',12, ...
         'markerfacecolor','b','markeredgecolor','b');    
    ha=xlabel('T (K)');
    hb=ylabel('D (m2/s)');
    set([ha hb],'fontsize',30,'fontWeight','bold');
    set(gca,'fontsize',30,'fontWeight','bold');
    NumTicks = 6;
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks));
    box on;
    savefig(fig3,'D2')

end



