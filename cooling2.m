clear;
clc;

%% paramters
s = 0.05;   %saturation paramter
Gamma = 2*pi*6.067*10^6;          %natural linewidth
Delta = -16*10^6;                 %detuning
k = 2*pi/(780*10^(-9));
m = 1.410*10^(-25);   %Rb-85
kb = 1.38*10^(-23);
T0 = 500*10^(-6);     %initial temperature
v0 = sqrt(3*kb*T0/m);
hbar = 6.63*10^(-34)/(2*pi);
N = 5000;              %total number of random walk steps
p = hbar*k/m;
time_axis = zeros(1,N);
temperature = zeros(1,N);

vx0 = v0/sqrt(3); 
vy0 = v0/sqrt(3);
vz0 = v0/sqrt(3);

%% main
% ax =@(t,vx) -hbar*k/m*Gamma/2*(s/(1+s+4*(Delta + k*vx)^2/Gamma^2) - s/(1+s+4*(Delta - k*vx)^2/Gamma^2));
% ay =@(t,vy) -hbar*k/m*Gamma/2*(s/(1+s+4*(Delta + k*vy)^2/Gamma^2) - s/(1+s+4*(Delta - k*vy)^2/Gamma^2));
% az =@(t,vz) -hbar*k/m*Gamma/2*(s/(1+s+4*(Delta + k*vz)^2/Gamma^2) - s/(1+s+4*(Delta - k*vz)^2/Gamma^2));
% 
% [tx,vx] = ode45(ax,tspan,vx0);
% [ty,vy] = ode45(ay,tspan,vy0);
% [tz,vz] = ode45(az,tspan,vz0);
% figure(1)
% hold on;
% T = m*(vx.^2 + vy.^2 + vz.^2)/(3*kb);
% plot(tx,T);
% hold off;
gamma_p = @(v) Gamma/2*s./(1+s+4*(Delta+k*v).^2/Gamma^2);

vx = vx0;
vy = vy0;
vz = vz0;
for i = 1:N-1
    v = [vx,-vx,vy,-vy,vz,-vz];
    temperature(i) = m*(vx^2 + vy^2 + vz^2)/(3*kb);
    rate = gamma_p(v);
    prob_norm = rate/sum(rate);
    random_num = rand();
    if  0 <= random_num && random_num < prob_norm(1)
        vx = vx - p;
        time_axis(i + 1) = time_axis(i) + 1/rate(1);
    elseif prob_norm(1) <= random_num && random_num < (prob_norm(1) + prob_norm(2))
        vx = vx + p;
        time_axis(i + 1) = time_axis(i) + 1/rate(2);
    elseif sum(prob_norm(1:2)) <= random_num && random_num < sum(prob_norm(1:3))
        vy = vy - p;
        time_axis(i + 1) = time_axis(i) + 1/rate(3);
    elseif sum(prob_norm(1:3)) <= random_num && random_num < sum(prob_norm(1:4))
        vy = vy + p;
        time_axis(i + 1) = time_axis(i) + 1/rate(4);
    elseif sum(prob_norm(1:4)) <= random_num && random_num < sum(prob_norm(1:5))
        vz = vz - p;
        time_axis(i + 1) = time_axis(i) + 1/rate(5);
    else 
        vz = vz + p;
        time_axis(i + 1) = time_axis(i) + 1/rate(6);
    end
end
temperature(N) = m*(vx^2 + vy^2 + vz^2)/(3*kb);

plot(time_axis*1000,temperature);
    