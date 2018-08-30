clear;
clc;

alpha0 = 0.5;
N = 5000;
T = 0:N;
x0 = 0;
x = zeros(1,N+1);
x(1) = x0;

for i=1:N
    alpha = rand();
    if alpha >= alpha0;
        x(i+1) = x(i) + 1;
    else
        x(i+1) = x(i) - 1;
    end
end

plot(T,x,'-o');