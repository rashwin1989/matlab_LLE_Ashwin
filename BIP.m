clear all;
close all;

T=[400:1:800];

dsp=[1 19;4 19;7 19;10 19;13 19;16 19;18 19];
[nPairs,nComp]=size(dsp);

d=load('petro.dat+');

kij_tmp=zeros(2,2);

for j=1:nPairs
    fprintf('%d\n',j);
    sp=dsp(j,:);
    Tc=d(sp,4);
    Pc=d(sp,5)*1e5;
    w =d(sp,6);    
    n=numel(w);
    tk=-sp';

    for i=1:numel(T)
        
        fprintf('%8.2f\n',T(i));
        [kij_tmp]=Matlab_BIP(T(i),Pc,Tc,w,tk);

        k(i,j)=kij_tmp(1,2);
    end
end

m=numel(T);

fig1 = figure(1);
hold on;      
plot(T,k(:,1),'-c','linewidth',4);
plot(T(1:10:m),k(1:10:m,1),'^c','markersize',12, ...
     'markerfacecolor','c','markeredgecolor','c');
plot(T,k(:,2),'-g','linewidth',4);
plot(T(1:10:m),k(1:10:m,2),'>g','markersize',12, ...
     'markerfacecolor','g','markeredgecolor','g');
plot(T,k(:,3),'-y','linewidth',4);
plot(T(1:10:m),k(1:10:m,3),'vy','markersize',12, ...
     'markerfacecolor','y','markeredgecolor','y');
plot(T,k(:,4),'-r','linewidth',4);
plot(T(1:10:m),k(1:10:m,4),'<r','markersize',12, ...
     'markerfacecolor','r','markeredgecolor','r');
plot(T,k(:,5),'-m','linewidth',4);
plot(T(1:10:m),k(1:10:m,5),'sm','markersize',12, ...
     'markerfacecolor','m','markeredgecolor','m');
plot(T,k(:,6),'-','Color',[0.8 0.5 0],'linewidth',4);
plot(T(1:10:m),k(1:10:m,6),'d','Color',[0.8 0.5 0],'markersize',12, ...
     'markerfacecolor',[0.8 0.5 0],'markeredgecolor',[0.8 0.5 0]);
plot(T,k(:,7),'-k','linewidth',4);
plot(T(1:10:m),k(1:10:m,7),'ok','markersize',12, ...
     'markerfacecolor','k','markeredgecolor','k');
hb=ylabel('k_{ij}');
ha=xlabel('T (K)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
NumTicks = 6;
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks));
box on;
savefig(fig1,'kij_all')