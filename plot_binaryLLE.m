clear all;
close all;

blle=load('binaryLLE.dat');
tT=blle(:,1);
single_ph=blle(:,2);
tx12=blle(:,3);
tx22=blle(:,4);

plotPt = 1;
for i=1:numel(tT)
    if(plotPt > 0)
        T(i,1) = tT(i,1);
        x12(i,1) = tx12(i,1);
        x22(i,1) = tx22(i,1);
    end
    if(single_ph(i) > 0)
        plotPt = -1;
    end
end

m=numel(T);

fig1 = figure(1);
hold on;      
plot(x12,T,'-b','linewidth',4);
plot(x12(1:m/10:m),T(1:m/10:m),'ob','markersize',12, ...
     'markerfacecolor','b','markeredgecolor','b');
plot(x22,T,'--b','linewidth',4);
plot(x22(1:m/10:m),T(1:m/10:m),'ob','markersize',12, ...
     'markerfacecolor','none','markeredgecolor','b');
ha=xlabel('X');
hb=ylabel('T (K)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
NumTicks = 6;
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks));
box on;
savefig(fig1,'x2_blle')