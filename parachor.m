clear all;
close all;

sp=[3 8 13 17 19];

P=30*1e6;

T=[600:1:800];

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
Vc=d(sp,8);    %cm3/mol
Vc=Vc/1000;
n=numel(w);
tk=-sp';

a0 = -72.765582; a1 = -9166.37732; a2 = -1.07786869;
a3 = 2197.28955; a4 = -111.984168; a5 = 20.4550431;
a6 = 55.9694357;
Pch=zeros(n,1);
for i=1:n-1
    % H = Vc(i)^(5/6)*Tc(i)^(1/4);
    % Pch(i) = a0 + a1*Vc(i) + a2*Tc(i) + a3*H + a4*H^2 + a5*H^3 +
    % a6/H; 
    Pch(i) = 69.9 + 2.3*MW(i);
end
%Pch(1)=2000;
sigma_water = 0.07*1e5/100;    %dynes/cm
Pch(n)=(sigma_water^0.25)*MW(n)/(1 - 0.001);

%d=load('Diaz.dat');
y0=load('y0.dat+');
y01=y0(1,:);
y02=y0(2,:);
x01=y01'.*invMW/(y01*invMW);
x02=y02'.*invMW/(y02*invMW);

x1=x01;
x2=x02;

for m12=[ 1 ]

    m1i=1/(1 + 1/m12);
    m2i=1-m1i; 
    c0=m1i/(x1'*MW)/(m1i/(x1'*MW) + m2i/(x2'*MW));

    display([num2str(c0) ': mixing ratio (wt) of crude vs. water is ' num2str(m1i) ]);
    
    for i=1:numel(T)
                
        [x1 x2 c]=Matlab_mLLE_new(c0,P,T(i),Pc,Tc,w,tk,x01,x02);

        [v1]=Matlab_volume(P,T(i),Pc,Tc,w,tk,x1);
        MW1 = x1'*MW;
        rho1 = 1e-6*MW1/v1;   %g/cm3
        [v2]=Matlab_volume(P,T(i),Pc,Tc,w,tk,x2);
        MW2 = x2'*MW;
        rho2 = 1e-6*MW2/v2;   %g/cm3

        sig=0;
        for j=1:n
            sig = sig + Pch(j)*(x1(j,1)*rho1/MW1 - x2(j,1)*rho2/MW2);
        end        
        sigma(i,1)=(sig^4)*1e-5*100;    %N/m from dynes/cm        
    end

    m=numel(T);

    fig1 = figure(1);
    hold on;
    plot(T,sigma(:,1),'-r','linewidth',4);
    plot(T(1:10:m),sigma(1:10:m,1),'or','markersize',12, ...
         'markerfacecolor','r','markeredgecolor','r');
    ha=xlabel('T (K)');
    hb=ylabel('\sigma (N/m)');
    set([ha hb],'fontsize',30,'fontWeight','bold');
    set(gca,'fontsize',30,'fontWeight','bold');
    NumTicks = 6;
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks));
    box on;
    savefig(fig1,'sigma')

end

