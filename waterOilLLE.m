clear all;
close all;

sp=[2 8 18 18 19];

%specific gravity of Shanghai oil: 0.9724
%20g oil in batch reactor, 20ml volume of oil, 80ml 8.5MPa N2
%charged = 7.6g N2; oil:N2 = 2.6

d=load('petro.dat+');
Tc=d(sp,4);
Pc=d(sp,5)*1e5;
w =d(sp,6);
MW=d(sp,7);
for i=1:numel(MW)
 invMW(i,1) = 1/MW(i,1);    
end

n=numel(w);

kij=zeros(1,n-1)*0 - 10000;

tk=-sp';

%d=load('Diaz.dat');
y0=load('y0.dat+');
y01=y0(1,:);
y02=y0(2,:);
x01=y01'.*invMW/(y01*invMW);
x02=y02'.*invMW/(y02*invMW);

x1=x01;
x2=x02;

for m12=[ 0.5 ]
    %for c0=[ 0.01:0.02:0.11 ]
    %for c0=[ 0.01:0.01:0.80 ]

    m1i=1/(1 + 1/m12);
    m2i=1-m1i; 
    c0=m1i/(x1'*MW)/(m1i/(x1'*MW) + m2i/(x2'*MW));

    display([num2str(c0) ': mixing ratio (wt) of crude vs. water is ' num2str(m1i) ]);

    P=30*1e6;

    T=[600:1:800];
    for i=1:numel(T)
                
        [x1 x2 c]=Matlab_mLLE_new(c0,P,T(i),Pc,Tc,w,tk,x01,x02);

        [v1]=Matlab_volume(P,T(i),Pc,Tc,w,tk,x1);
        [v2]=Matlab_volume(P,T(i),Pc,Tc,w,tk,x2);

        rho1 = 1e-3*(x1'*MW)/v1;
        rho2 = 1e-3*(x2'*MW)/v2;

        C1(i)=c;
        m1(i)=c*(x1'*MW)/(c*(x1'*MW) + (1 - c)*(x2'*MW));
        y1(:,i)=x1.*MW/(x1'*MW);
        y2(:,i)=x2.*MW/(x2'*MW);
        He(:,i)=(rho1/rho2)*(y1(:,i)./y2(:,i));
        X1(:,i)=x1;
        X2(:,i)=x2;
    end

    m=numel(T);

    switch n

      case 2

        fig1 = figure(1);
        hold on;      
        plot(y1(1,:),T,'-k','linewidth',4);
        plot(y1(1,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','k','markeredgecolor','k');
        plot(y2(1,:),T,'--k','linewidth',4);
        plot(y2(1,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','k');
        plot(y1(2,:),T,'-b','linewidth',4);
        plot(y1(2,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','b','markeredgecolor','b');
        plot(y2(2,:),T,'--b','linewidth',4);
        plot(y2(2,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','b');
        ha=xlabel('Y');
        hb=ylabel('T (K)');
        set([ha hb],'fontsize',30,'fontWeight','bold');
        set(gca,'fontsize',30,'fontWeight','bold');
        NumTicks = 6;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        box on;
        savefig(fig1,'y_all')

      case 3

        fig1 = figure(1);
        hold on;      
        plot(y1(1,:),T,'-r','linewidth',4);
        plot(y1(1,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','r','markeredgecolor','r');
        plot(y2(1,:),T,'--r','linewidth',4);
        plot(y2(1,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','r');
        plot(y1(2,:),T,'-k','linewidth',4);
        plot(y1(2,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','k','markeredgecolor','k');
        plot(y2(2,:),T,'--k','linewidth',4);
        plot(y2(2,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','k');
        plot(y1(3,:),T,'-b','linewidth',4);
        plot(y1(3,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','b','markeredgecolor','b');
        plot(y2(3,:),T,'--b','linewidth',4);
        plot(y2(3,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','b');
        ha=xlabel('Y');
        hb=ylabel('T (K)');
        set([ha hb],'fontsize',30,'fontWeight','bold');
        set(gca,'fontsize',30,'fontWeight','bold');
        NumTicks = 6;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        box on;
        savefig(fig1,'y_all')

      case 4

        fig1 = figure(1);
        hold on;
        plot(y1(1,:),T,'-g','linewidth',4);
        plot(y1(1,1:10:m),T(1:10:m),'^g','markersize',12, ...
             'markerfacecolor','g','markeredgecolor','g');
        plot(y2(1,:),T,'--g','linewidth',4);
        plot(y2(1,1:10:m),T(1:10:m),'^g','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','g');
        plot(y1(2,:),T,'-r','linewidth',4);
        plot(y1(2,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','r','markeredgecolor','r');
        plot(y2(2,:),T,'--r','linewidth',4);
        plot(y2(2,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','r');
        plot(y1(3,:),T,'-k','linewidth',4);
        plot(y1(3,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','k','markeredgecolor','k');
        plot(y2(3,:),T,'--k','linewidth',4);
        plot(y2(3,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','k');
        plot(y1(4,:),T,'-b','linewidth',4);
        plot(y1(4,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','b','markeredgecolor','b');
        plot(y2(4,:),T,'--b','linewidth',4);
        plot(y2(4,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','b');
        ha=xlabel('Y');
        hb=ylabel('T (K)');
        set([ha hb],'fontsize',30,'fontWeight','bold');
        set(gca,'fontsize',30,'fontWeight','bold');
        NumTicks = 6;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        box on;
        savefig(fig1,'y_all')

      case 5

        fig1 = figure(1);
        hold on;
        plot(y1(1,:),T,'-g','linewidth',4);
        plot(y1(1,1:10:m),T(1:10:m),'^g','markersize',12, ...
             'markerfacecolor','g','markeredgecolor','g');
        plot(y2(1,:),T,'--g','linewidth',4);
        plot(y2(1,1:10:m),T(1:10:m),'^g','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','g');
        plot(y1(2,:),T,'-y','linewidth',4);
        plot(y1(2,1:10:m),T(1:10:m),'>y','markersize',12, ...
             'markerfacecolor','y','markeredgecolor','y');
        plot(y2(2,:),T,'--y','linewidth',4);
        plot(y2(2,1:10:m),T(1:10:m),'>y','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','y');
        plot(y1(3,:),T,'-r','linewidth',4);
        plot(y1(3,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','r','markeredgecolor','r');
        plot(y2(3,:),T,'--r','linewidth',4);
        plot(y2(3,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','r');        
        plot(y1(4,:),T,'-k','linewidth',4);
        plot(y1(4,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','k','markeredgecolor','k');
        plot(y2(4,:),T,'--k','linewidth',4);
        plot(y2(4,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','k');
        plot(y1(5,:),T,'-b','linewidth',4);
        plot(y1(5,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','b','markeredgecolor','b');
        plot(y2(5,:),T,'--b','linewidth',4);
        plot(y2(5,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','b');
        ha=xlabel('Y');
        hb=ylabel('T (K)');
        set([ha hb],'fontsize',30,'fontWeight','bold');
        set(gca,'fontsize',30,'fontWeight','bold');
        NumTicks = 6;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        box on;
        savefig(fig1,'y_all')        

      case 6

        fig1 = figure(1);
        hold on;
        plot(y1(1,:),T,'-g','linewidth',4);
        plot(y1(1,1:10:m),T(1:10:m),'^g','markersize',12, ...
             'markerfacecolor','g','markeredgecolor','g');
        plot(y2(1,:),T,'--g','linewidth',4);
        plot(y2(1,1:10:m),T(1:10:m),'^g','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','g');
        plot(y1(2,:),T,'-y','linewidth',4);
        plot(y1(2,1:10:m),T(1:10:m),'>y','markersize',12, ...
             'markerfacecolor','y','markeredgecolor','y');
        plot(y2(2,:),T,'--y','linewidth',4);
        plot(y2(2,1:10:m),T(1:10:m),'>y','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','y');
        plot(y1(3,:),T,'-r','linewidth',4);
        plot(y1(3,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','r','markeredgecolor','r');
        plot(y2(3,:),T,'--r','linewidth',4);
        plot(y2(3,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','r');
        plot(y1(4,:),T,'-m','linewidth',4);
        plot(y1(4,1:10:m),T(1:10:m),'dm','markersize',12, ...
             'markerfacecolor','m','markeredgecolor','m');
        plot(y2(4,:),T,'--m','linewidth',4);
        plot(y2(4,1:10:m),T(1:10:m),'dm','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','m');
        plot(y1(5,:),T,'-k','linewidth',4);
        plot(y1(5,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','k','markeredgecolor','k');
        plot(y2(5,:),T,'--k','linewidth',4);
        plot(y2(5,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','k');
        plot(y1(6,:),T,'-b','linewidth',4);
        plot(y1(6,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','b','markeredgecolor','b');
        plot(y2(6,:),T,'--b','linewidth',4);
        plot(y2(6,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','b');
        ha=xlabel('Y');
        hb=ylabel('T (K)');
        set([ha hb],'fontsize',30,'fontWeight','bold');
        set(gca,'fontsize',30,'fontWeight','bold');
        NumTicks = 6;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        box on;
        savefig(fig1,'y_all')

        case 10

        fig1 = figure(1);
        hold on;
        plot(y1(1,:),T,'-g','linewidth',4);
        plot(y1(1,1:10:m),T(1:10:m),'^g','markersize',12, ...
             'markerfacecolor','g','markeredgecolor','g');
        plot(y2(1,:),T,'--g','linewidth',4);
        plot(y2(1,1:10:m),T(1:10:m),'^g','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','g');
        plot(y1(4,:),T,'-y','linewidth',4);
        plot(y1(4,1:10:m),T(1:10:m),'>y','markersize',12, ...
             'markerfacecolor','y','markeredgecolor','y');
        plot(y2(4,:),T,'--y','linewidth',4);
        plot(y2(4,1:10:m),T(1:10:m),'>y','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','y');
        plot(y1(6,:),T,'-r','linewidth',4);
        plot(y1(6,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','r','markeredgecolor','r');
        plot(y2(6,:),T,'--r','linewidth',4);
        plot(y2(6,1:10:m),T(1:10:m),'vr','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','r');
        plot(y1(8,:),T,'-m','linewidth',4);
        plot(y1(8,1:10:m),T(1:10:m),'dm','markersize',12, ...
             'markerfacecolor','m','markeredgecolor','m');
        plot(y2(8,:),T,'--m','linewidth',4);
        plot(y2(8,1:10:m),T(1:10:m),'dm','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','m');
        plot(y1(9,:),T,'-k','linewidth',4);
        plot(y1(9,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','k','markeredgecolor','k');
        plot(y2(9,:),T,'--k','linewidth',4);
        plot(y2(9,1:10:m),T(1:10:m),'sk','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','k');
        plot(y1(10,:),T,'-b','linewidth',4);
        plot(y1(10,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','b','markeredgecolor','b');
        plot(y2(10,:),T,'--b','linewidth',4);
        plot(y2(10,1:10:m),T(1:10:m),'ob','markersize',12, ...
             'markerfacecolor','none','markeredgecolor','b');
        ha=xlabel('Y');
        hb=ylabel('T (K)');
        set([ha hb],'fontsize',30,'fontWeight','bold');
        set(gca,'fontsize',30,'fontWeight','bold');
        NumTicks = 6;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        box on;
        savefig(fig1,'y_all')

    end


    fig2 = figure(2);
    hold on;
    plot(T,m1,'-k','linewidth',4);
    plot(T(1:10:m),m1(1:10:m),'ok','markersize',12, ...
         'markerfacecolor','k','markeredgecolor','k');
    ha=xlabel('T (K)');
    hb=ylabel('Yoil');
    set([ha hb],'fontsize',30,'fontWeight','bold');
    set(gca,'fontsize',30,'fontWeight','bold');
    NumTicks = 6;
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks));
    box on;
    savefig(fig2,'Yoil')

    fig3 = figure(3);
    hold on;
    plot(T,He(1,:),'-g','linewidth',4);
    plot(T(1:10:m),He(1,1:10:m),'og','markersize',12, ...
         'markerfacecolor','g','markeredgecolor','g');
    ha=xlabel('T (K)');
    hb=ylabel('He');
    set([ha hb],'fontsize',30,'fontWeight','bold');
    set(gca,'fontsize',30,'fontWeight','bold');
    NumTicks = 6;
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks));
    box on;
    savefig(fig3,'He_lightOil')

end

