clear all;
close all;

sp=[3 8 13 17 19];

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
y0=load('y0.dat+');
y01=y0(1,:);
y02=y0(2,:);
x01=y01'.*invMW/(y01*invMW);
x02=y02'.*invMW/(y02*invMW);

x1=x01;
x2=x02;

Dij1_tmp=zeros(n,n);
Dij2_tmp=zeros(n,n);

for m12=[ 1 ]
    
    m1i=1/(1 + 1/m12);
    m2i=1-m1i; 
    c0=m1i/(x1'*MW)/(m1i/(x1'*MW) + m2i/(x2'*MW));

    display([num2str(c0) ': mixing ratio (wt) of crude vs. water is ' num2str(m1i) ]);

    P=30*1e6;

    T=[600:1:760];
    for i=1:numel(T)
                
        [x1 x2 c]=Matlab_mLLE_new(c0,P,T(i),Pc,Tc,w,tk,x01,x02);

        [v1]=Matlab_volume(P,T(i),Pc,Tc,w,tk,x1);
        [v2]=Matlab_volume(P,T(i),Pc,Tc,w,tk,x2);

        [mu1_tmp lambda1_tmp]=Matlab_vis_cond(P,T(i),Pc,Tc,Vc,w, ...
                                               MW,kappa,dm,x1,v1,Tb,SG,H8);
        [mu2_tmp lambda2_tmp]=Matlab_vis_cond(P,T(i),Pc,Tc,Vc,w, ...
                                               MW,kappa,dm,x2,v2, ...
                                               Tb,SG,H8);

        [Dij1_tmp]=Matlab_Dij_tlsm_wk(P,T(i),Pc,Tc,Vc,w,tk,MW,x1);
        [Dij2_tmp]=Matlab_Dij_tlsm_wk(P,T(i),Pc,Tc,Vc,w,tk,MW,x2);

        for k=1:n
            Di1_tmp(k) = 0;
            Di2_tmp(k) = 0;
            for l=1:n
                if (k ~= l)
                    Di1_tmp(k) = Di1_tmp(k) + x1(l)/Dij1_tmp(k,l);
                    Di2_tmp(k) = Di2_tmp(k) + x2(l)/Dij2_tmp(k,l);
                end
            end
            Di1_tmp(k) = (1 - x1(k))/Di1_tmp(k);
            Di2_tmp(k) = (1 - x2(k))/Di2_tmp(k);
        end

        rho1(i,1) = 1e-3*(x1'*MW)/v1;
        mu1(i,1) = mu1_tmp;
        lambda1(i,1) = lambda1_tmp;
        Di1(i,:) = Di1_tmp;
        
        rho2(i,1) = 1e-3*(x2'*MW)/v2;
        mu2(i,1) = mu2_tmp;
        lambda2(i,1) = lambda2_tmp;
        Di2(i,:) = Di2_tmp;
                                               
    end

    m=numel(T);

    switch n      

      case 5

        fig1 = figure(1);
        hold on;
        plot(T,rho1(:,1),'-r','linewidth',4);
        plot(T(1:10:m),rho1(1:10:m,1),'sr','markersize',12, ...
             'markerfacecolor','r','markeredgecolor','r');
        plot(T,rho2(:,1),'-b','linewidth',4);
        plot(T(1:10:m),rho2(1:10:m,1),'ob','markersize',12, ...
             'markerfacecolor','b','markeredgecolor','b');
        ha=xlabel('T (K)');
        hb=ylabel('\rho (kg/m3)');
        set([ha hb],'fontsize',30,'fontWeight','bold');
        set(gca,'fontsize',30,'fontWeight','bold');
        NumTicks = 6;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        box on;
        savefig(fig1,'rho')

        % fig2 = figure(2);
        % hold on;
        % plot(T(1:100),mu1(1:100,1),'-r','linewidth',4);
        % plot(T(1:10:100),mu1(1:10:100,1),'sr','markersize',12, ...
        %      'markerfacecolor','r','markeredgecolor','r');
        % plot(T(1:100),mu2(1:100,1),'-b','linewidth',4);
        % plot(T(1:10:100),mu2(1:10:100,1),'ob','markersize',12, ...
        %      'markerfacecolor','b','markeredgecolor','b');
        % ha=xlabel('T (K)');
        % hb=ylabel('\mu (Pa-s)');
        % set([ha hb],'fontsize',30,'fontWeight','bold');
        % set(gca,'fontsize',30,'fontWeight','bold');
        % NumTicks = 6;
        % L = get(gca,'YLim');
        % set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        % box on;
        % savefig(fig2,'mu')

        fig3 = figure(3);
        hold on;
        plot(T,Di1(:,1),'-g','linewidth',4);
        plot(T(1:10:m),Di1(1:10:m,1),'^g','markersize',12, ...
             'markerfacecolor','g','markeredgecolor','g');
        plot(T,Di2(:,1),'--g','linewidth',4);
        plot(T(1:10:m),Di2(1:10:m,1),'^g','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','g');
        plot(T,Di1(:,3),'-r','linewidth',4);
        plot(T(1:10:m),Di1(1:10:m,3),'sr','markersize',12, ...
             'markerfacecolor','r','markeredgecolor','r');
        plot(T,Di2(:,3),'--r','linewidth',4);
        plot(T(1:10:m),Di2(1:10:m,3),'sr','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','r');
        plot(T,Di1(:,5),'-b','linewidth',4);
        plot(T(1:10:m),Di1(1:10:m,5),'ob','markersize',12, ...
             'markerfacecolor','b','markeredgecolor','b');
        plot(T,Di2(:,5),'--b','linewidth',4);
        plot(T(1:10:m),Di2(1:10:m,5),'ob','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','b');
        ha=xlabel('T (K)');
        hb=ylabel('D (m2/s)');
        set([ha hb],'fontsize',30,'fontWeight','bold');
        set(gca,'fontsize',30,'fontWeight','bold');
        NumTicks = 6;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        box on;
        savefig(fig3,'D')

    end

end

