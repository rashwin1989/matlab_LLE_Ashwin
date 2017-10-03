function [A,B] = calculate_A_B(w)

N = numel(w);
A = zeros(N,1);
B = zeros(N,1);

for i=1:N
    A(i) = -1e-16/(1.2326 + 1.3757*w(i));
    B(i) = 1e-16/(0.9051 + 1.5410*w(i));
end