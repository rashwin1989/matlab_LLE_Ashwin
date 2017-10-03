clear all;
close all;

b1 = 23.55379;
b2 = -3.69839;
b3 = 0.005766;

P=[20:5:30];
T=[300:5:800];
for j=1:numel(P)
    for i=1:numel(T)
        mu(i,j) = 1e-3*exp(b1 + b2*log(T(i)) + b3*P(j));
    end
end

m = numel(T);

fig1 = figure(1);
hold on;      
plot(T,mu(:,1),'-g','linewidth',4);
plot(T(1:10:m),mu(1:10:m,1),'sg','markersize',12, ...
     'markerfacecolor','g','markeredgecolor','g');
plot(T,mu(:,2),'-b','linewidth',4);
plot(T(1:10:m),mu(1:10:m,2),'db','markersize',12, ...
     'markerfacecolor','b','markeredgecolor','b');
plot(T,mu(:,3),'-r','linewidth',4);
plot(T(1:10:m),mu(1:10:m,3),'or','markersize',12, ...
     'markerfacecolor','r','markeredgecolor','r');
ha=xlabel('T (K)');
hb=ylabel('mu (Pa-s)');
set([ha hb],'fontsize',30,'fontWeight','bold');
set(gca,'fontsize',30,'fontWeight','bold');
NumTicks = 6;
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks));
box on;
savefig(fig1,'mu_bitumen')