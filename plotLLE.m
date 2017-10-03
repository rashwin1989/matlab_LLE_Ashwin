clear all;
close all;

lle=load('LLE.dat');
T=lle(:,1);
m1=lle(:,2);
single_ph=lle(:,3);

y1=load('y1.dat');
y2=load('y2.dat');

[m,n]=size(y1);

switch n

    case 2

      fig1 = figure(1);
      hold on;      
      plot(y1(:,1),T,'-k','linewidth',4);
      plot(y1(1:10:m,1),T(1:10:m),'sk','markersize',12, ...
           'markerfacecolor','k','markeredgecolor','k');
      plot(y2(:,1),T,'--k','linewidth',4);
      plot(y2(1:10:m,1),T(1:10:m),'sk','markersize',12, ...
           'markerfacecolor','none','markeredgecolor','k');
      plot(y1(:,2),T,'-b','linewidth',4);
      plot(y1(1:10:m,2),T(1:10:m),'ob','markersize',12, ...
           'markerfacecolor','b','markeredgecolor','b');
      plot(y2(:,2),T,'--b','linewidth',4);
      plot(y2(1:10:m,2),T(1:10:m),'ob','markersize',12, ...
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
      plot(y1(:,1),T,'-r','linewidth',4);
      plot(y1(1:10:m,1),T(1:10:m),'vr','markersize',12, ...
           'markerfacecolor','r','markeredgecolor','r');
      plot(y2(:,1),T,'--r','linewidth',4);
      plot(y2(1:10:m,1),T(1:10:m),'vr','markersize',12, ...
           'markerfacecolor','none','markeredgecolor','r');
      plot(y1(:,2),T,'-k','linewidth',4);
      plot(y1(1:10:m,2),T(1:10:m),'sk','markersize',12, ...
           'markerfacecolor','k','markeredgecolor','k');
      plot(y2(:,2),T,'--k','linewidth',4);
      plot(y2(1:10:m,2),T(1:10:m),'sk','markersize',12, ...
           'markerfacecolor','none','markeredgecolor','k');
      plot(y1(:,3),T,'-b','linewidth',4);
      plot(y1(1:10:m,3),T(1:10:m),'ob','markersize',12, ...
           'markerfacecolor','b','markeredgecolor','b');
      plot(y2(:,3),T,'--b','linewidth',4);
      plot(y2(1:10:m,3),T(1:10:m),'ob','markersize',12, ...
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
      plot(y1(:,1),T,'-g','linewidth',4);
      plot(y1(1:10:m,1),T(1:10:m),'^g','markersize',12, ...
           'markerfacecolor','g','markeredgecolor','g');
      plot(y2(:,1),T,'--g','linewidth',4);
      plot(y2(1:10:m,1),T(1:10:m),'^g','markersize',12, ...
           'markerfacecolor','none','markeredgecolor','g');
      plot(y1(:,2),T,'-r','linewidth',4);
      plot(y1(1:10:m,2),T(1:10:m),'vr','markersize',12, ...
           'markerfacecolor','r','markeredgecolor','r');
      plot(y2(:,2),T,'--r','linewidth',4);
      plot(y2(1:10:m,2),T(1:10:m),'vr','markersize',12, ...
           'markerfacecolor','none','markeredgecolor','r');
      plot(y1(:,3),T,'-k','linewidth',4);
      plot(y1(1:10:m,3),T(1:10:m),'sk','markersize',12, ...
           'markerfacecolor','k','markeredgecolor','k');
      plot(y2(:,3),T,'--k','linewidth',4);
      plot(y2(1:10:m,3),T(1:10:m),'sk','markersize',12, ...
           'markerfacecolor','none','markeredgecolor','k');
      plot(y1(:,4),T,'-b','linewidth',4);
      plot(y1(1:10:m,4),T(1:10:m),'ob','markersize',12, ...
           'markerfacecolor','b','markeredgecolor','b');
      plot(y2(:,4),T,'--b','linewidth',4);
      plot(y2(1:10:m,4),T(1:10:m),'ob','markersize',12, ...
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