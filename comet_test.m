clear;
clc;

U0 = 1;
w0 = 1.8 * 10^(-6);   %waist of trap, unit: m
lambda_trap = 828 * 10^(-9); %wavelength of trap laser, unit:m
z_R = pi*w0^2/lambda_trap;      %Rayleigh length, unit:m
U =@(rvec) -U0/(1 + (rvec(3)/z_R)^2) * exp(-2*(rvec(1)^2 + rvec(2)^2)/(w0^2*(1 + rvec(3)^2/z_R^2)));

N = 1000;
x = transpose(linspace(-w0,w0,N));
y = transpose(linspace(-w0,w0,N));
z = transpose(linspace(-z_R,z_R,N));
Coordinate = [x,y,z];
Potential = zeros(N,1);

for i = 1:N
    Potential(i) = U(Coordinate(i,:));
end
surf(x,y,z,Potential);
