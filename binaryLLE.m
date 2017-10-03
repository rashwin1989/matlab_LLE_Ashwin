clear all;
close all;

P=30*1e6;
T=[600:1:800];

x1=[1.0;0];
x2=[0;1.0];

dsp=[1 19;4 19;7 19;10 19;13 19;16 19;18 19];
[nPairs,nComp]=size(dsp);

d=load('petro.dat+');

for j=1:nPairs
    fprintf('%d\n',j);
    sp=dsp(j,:);
    Tc=d(sp,4);
    Pc=d(sp,5)*1e5;
    w =d(sp,6);
    MW=d(sp,7);
    for i=1:numel(MW)
        invMW(i,1) = 1/MW(i,1);    
    end
    n=numel(w);
    tk=-sp';

    for i=1:numel(T)
        
        fprintf('%8.2f\n',T(i));
        [x1 x2 n_miscible]=Matlab_binaryLLE(P,T(i),Pc,Tc,w,tk);

        if (abs(x1(1,1)-x2(1,1)) < 1e-6)
            x1(1,1)=NaN;x1(2,1)=NaN;x2(1,1)=NaN;x2(2,1)=NaN;
            X1(:,i,j)=x1;
            X2(:,i,j)=x2;
            y1(:,i,j)=x1;
            y2(:,i,j)=x2;
        else
            y1(:,i,j)=x1.*MW/(x1'*MW);
            y2(:,i,j)=x2.*MW/(x2'*MW);
            X1(:,i,j)=x1;
            X2(:,i,j)=x2;
        end
    end
end

m=numel(T);

fig1 = figure(1);
hold on;      
plot(y1(1,:,1),T,'-c','linewidth',4);
plot(y1(1,1:10:m,1),T(1:10:m),'^c','markersize',12, ...
     'markerfacecolor','c','markeredgecolor','c');
plot(y2(1,:,1),T,'-c','linewidth',4);
plot(y2(1,1:10:m,1),T(1:10:m),'^c','markersize',12, ...
     'markerfacecolor','c','markeredgecolor','c');
plot(y1(1,:,2),T,'-g','linewidth',4);
plot(y1(1,1:10:m,2),T(1:10:m),'>g','markersize',12, ...
     'markerfacecolor','g','markeredgecolor','g');
plot(y2(1,:,2),T,'-g','linewidth',4);
plot(y2(1,1:10:m,2),T(1:10:m),'>g','markersize',12, ...
     'markerfacecolor','g','markeredgecolor','g');
plot(y1(1,:,3),T,'-y','linewidth',4);
plot(y1(1,1:10:m,3),T(1:10:m),'vy','markersize',12, ...
     'markerfacecolor','y','markeredgecolor','y');
plot(y2(1,:,3),T,'-y','linewidth',4);
plot(y2(1,1:10:m,3),T(1:10:m),'vy','markersize',12, ...
     'markerfacecolor','y','markeredgecolor','y');
plot(y1(1,:,4),T,'-r','linewidth',4);
plot(y1(1,1:10:m,4),T(1:10:m),'<r','markersize',12, ...
     'markerfacecolor','r','markeredgecolor','r');
plot(y2(1,:,4),T,'-r','linewidth',4);
plot(y2(1,1:10:m,4),T(1:10:m),'<r','markersize',12, ...
     'markerfacecolor','r','markeredgecolor','r');
plot(y1(1,:,5),T,'-m','linewidth',4);
plot(y1(1,1:10:m,5),T(1:10:m),'sm','markersize',12, ...
     'markerfacecolor','m','markeredgecolor','m');
plot(y2(1,:,5),T,'-m','linewidth',4);
plot(y2(1,1:10:m,5),T(1:10:m),'sm','markersize',12, ...
     'markerfacecolor','m','markeredgecolor','m');
plot(y1(1,:,6),T,'-','Color',[0.8 0.5 0],'linewidth',4);
plot(y1(1,1:10:m,6),T(1:10:m),'d','Color',[0.8 0.5 0],'markersize',12, ...
     'markerfacecolor',[0.8 0.5 0],'markeredgecolor',[0.8 0.5 0]);
plot(y2(1,:,6),T,'-','Color',[0.8 0.5 0],'linewidth',4);
plot(y2(1,1:10:m,6),T(1:10:m),'d','Color',[0.8 0.5 0],'markersize',12, ...
     'markerfacecolor',[0.8 0.5 0],'markeredgecolor',[0.8 0.5 0]);
plot(y1(1,:,7),T,'-k','linewidth',4);
plot(y1(1,1:10:m,7),T(1:10:m),'ok','markersize',12, ...
     'markerfacecolor','k','markeredgecolor','k');
plot(y2(1,:,7),T,'-k','linewidth',4);
plot(y2(1,1:10:m,7),T(1:10:m),'ok','markersize',12, ...
     'markerfacecolor','k','markeredgecolor','k');
ha=xlabel('Y');
hb=ylabel('T (K)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
NumTicks = 6;
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks));
box on;
savefig(fig1,'y_binary_all')


fig2 = figure(2);
hold on;      
plot(X1(1,:,1),T,'-c','linewidth',4);
plot(X1(1,1:10:m,1),T(1:10:m),'^c','markersize',12, ...
     'markerfacecolor','c','markeredgecolor','c');
plot(X2(1,:,1),T,'-c','linewidth',4);
plot(X2(1,1:10:m,1),T(1:10:m),'^c','markersize',12, ...
     'markerfacecolor','c','markeredgecolor','c');
plot(X1(1,:,2),T,'-g','linewidth',4);
plot(X1(1,1:10:m,2),T(1:10:m),'>g','markersize',12, ...
     'markerfacecolor','g','markeredgecolor','g');
plot(X2(1,:,2),T,'-g','linewidth',4);
plot(X2(1,1:10:m,2),T(1:10:m),'>g','markersize',12, ...
     'markerfacecolor','g','markeredgecolor','g');
plot(X1(1,:,3),T,'-y','linewidth',4);
plot(X1(1,1:10:m,3),T(1:10:m),'vy','markersize',12, ...
     'markerfacecolor','y','markeredgecolor','y');
plot(X2(1,:,3),T,'-y','linewidth',4);
plot(X2(1,1:10:m,3),T(1:10:m),'vy','markersize',12, ...
     'markerfacecolor','y','markeredgecolor','y');
plot(X1(1,:,4),T,'-r','linewidth',4);
plot(X1(1,1:10:m,4),T(1:10:m),'<r','markersize',12, ...
     'markerfacecolor','r','markeredgecolor','r');
plot(X2(1,:,4),T,'-r','linewidth',4);
plot(X2(1,1:10:m,4),T(1:10:m),'<r','markersize',12, ...
     'markerfacecolor','r','markeredgecolor','r');
plot(X1(1,:,5),T,'-m','linewidth',4);
plot(X1(1,1:10:m,5),T(1:10:m),'sm','markersize',12, ...
     'markerfacecolor','m','markeredgecolor','m');
plot(X2(1,:,5),T,'-m','linewidth',4);
plot(X2(1,1:10:m,5),T(1:10:m),'sm','markersize',12, ...
     'markerfacecolor','m','markeredgecolor','m');
plot(X1(1,:,6),T,'-','Color',[0.8 0.5 0],'linewidth',4);
plot(X1(1,1:10:m,6),T(1:10:m),'d','Color',[0.8 0.5 0],'markersize',12, ...
     'markerfacecolor',[0.8 0.5 0],'markeredgecolor',[0.8 0.5 0]);
plot(X2(1,:,6),T,'-','Color',[0.8 0.5 0],'linewidth',4);
plot(X2(1,1:10:m,6),T(1:10:m),'d','Color',[0.8 0.5 0],'markersize',12, ...
     'markerfacecolor',[0.8 0.5 0],'markeredgecolor',[0.8 0.5 0]);
plot(X1(1,:,7),T,'-k','linewidth',4);
plot(X1(1,1:10:m,7),T(1:10:m),'ok','markersize',12, ...
     'markerfacecolor','k','markeredgecolor','k');
plot(X2(1,:,7),T,'-k','linewidth',4);
plot(X2(1,1:10:m,7),T(1:10:m),'ok','markersize',12, ...
     'markerfacecolor','k','markeredgecolor','k');
ha=xlabel('X');
hb=ylabel('T (K)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
NumTicks = 6;
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks));
box on;
savefig(fig2,'x_binary_all')





