function [c] = calculate_ci(T,Pc,Tc,w,A,B)

N = numel(w);
c = zeros(N,1);
a = zeros(N,1);
b = zeros(N,1);
[a,b] = Matlab_a_b(T,Pc,Tc,w);

for i=1:N
    ti = 1.0 - T/Tc(i);
    c(i) = a(i)*b(i)^(2/3)*(A(i)*ti + B(i));
end
