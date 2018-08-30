function [rvec10, rvec20, v10, v20] = initial_r_v(w0,R_C,mbeta)
%Determine the initial position and velocity of atoms

r10 = rand()*w0/2;
theta10 = rand()*2*pi;
x10 = r10 * cos(theta10);
y10 = r10 * sin(theta10);
z10 = 0;                     %initial position of atom 1

r20 = rand()*w0/2;
theta20 = rand()*2*pi;
x20 = r20 * cos(theta20);
y20 = r20 * sin(theta20);
z20 = 0;                     %initial position of atom 2

R0 = sqrt((x10- x20)^2 + (y10 - y20)^2);
while R0 <= R_C || R0 >= w0
    r10 = rand()*w0/2;
    theta10 = rand()*2*pi;
    x10 = r10 * cos(theta10);
    y10 = r10 * sin(theta10);
    r20 = rand()*w0/2;
    theta20 = rand()*2*pi;
    x20 = r20 * cos(theta20);
    y20 = r20 * sin(theta20);
    R0 = sqrt((x10- x20)^2 + (y10 - y20)^2);
end          %guarantee the distance is larger than R_C and smaller than w0;

rvec10 = [x10,y10,z10];
rvec20 = [x20,y20,z20];

v1x0 = sqrt(1/mbeta)*randn();
v1y0 = sqrt(1/mbeta)*randn();
v1z0 = sqrt(1/mbeta)*randn();      %initial velocity of atom 2

v2x0 = sqrt(1/mbeta)*randn();
v2y0 = sqrt(1/mbeta)*randn();
v2z0 = sqrt(1/mbeta)*randn();       %initial velocity of atom 1

v10 = [v1x0,v1y0,v1z0];
v20 = [v2x0,v2y0,v2z0];

end

