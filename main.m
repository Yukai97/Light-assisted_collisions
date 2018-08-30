clear;
clc;

%% parameters
h = 6.63*10^(-34);       %Planck constant, unit: J*s
hbar = h/(2*pi);
g = 9.8;    %gravitational accelaration, unit:m/s^2
kb = 1.38*10^(-23);      %Boltzmann constant, unit: J/K

m = 1.410*10^(-25);   %Rb-85, unit: kg
mu = m/2;             %reduced mass
Gamma = 2*pi*6.067*10^6;          %natural linewidth, unit:Hz
k = 2*pi/(780*10^(-9));           %unit: m^(-1)
TD = hbar*Gamma/(2*kb);           %unit: K

T0 = 150*10^(-6);     %initial temperature, unit: K
mbeta = m/(kb * T0);        %m * beta in Maxwell distribution, unit: kg/J

s_cooling = 0.05;   %saturation paramter of cooling light
Delta_Cooling = -16*10^6;         %detuning of cooling light, unit: Hz
vp = hbar*k/m;          %recoil velocity, unit:m/s
gamma_p = @(v) Gamma/2*s_cooling./(1+s_cooling+4*(Delta_Cooling+k*v).^2/Gamma^2);   %scattering rate, unit:s^(-1)

Delta_Collision = -45*10^6;       %detuning of collision light, unit: Hz
s_collision = 0.1;
Rabi_collision_square = s_collision * Gamma^2/2;    %The square of Rabi frequency

U0 = hbar*85*10^6/(kb*TD);         %trap depth, unit:k_b * T_D
w0 = 1.8 * 10^(-6);   %waist of trap, unit: m
lambda_trap = 828 * 10^(-9); %wavelength of trap laser, unit:m
z_R = pi*w0^2/lambda_trap;      %Rayleigh length, unit:m
U =@(rvec) -U0/(1 + (rvec(3)/z_R)^2) * exp(-2*(rvec(1)^2 + rvec(2)^2)/(w0^2*(1 + rvec(3)^2/z_R^2)));     % trap function
Fx =@(rvec) -4*U0*exp(-2*(rvec(1)^2 + rvec(2)^2)/(w0^2*(1 + rvec(3)^2/z_R^2)))/(w0^2*(1 + (rvec(3)/z_R)^2))*rvec(1);     
Fy =@(rvec) -4*U0*exp(-2*(rvec(1)^2 + rvec(2)^2)/(w0^2*(1 + rvec(3)^2/z_R^2)))/(w0^2*(1 + (rvec(3)/z_R)^2))*rvec(2);
Fz =@(rvec) Trap_Fz(rvec,w0,U0,z_R);    
F = {Fx,Fy,Fz};                    %force corresponding to trap

C3 = -20.13*4.3597*(5.29177)^3/(kb*10^23*TD)*10^(-28);   %molecular potential constant, unit: kb*TD/m^3
Umol =@(rvec) -C3/norm(rvec)^3;
Fmolx = @(rvec) -3*C3/norm(rvec)^5*rvec(1);
Fmoly = @(rvec) -3*C3/norm(rvec)^5*rvec(2);
Fmolz = @(rvec) -3*C3/norm(rvec)^5*rvec(3);       
Fmol = {Fmolx,Fmoly,Fmolz};               %force corresponding to molecular potentials
R_Condon = (C3/(hbar*Delta_Collision)*kb*TD)^(1/3);     %Condon point
Reaction_D = 1/100 * R_Condon;
P_LZ = @(vr) exp(2*pi*(hbar*Rabi_collision_square/(kb*TD))*(R_Condon^4/C3)/(3*vr));

N = 50000;              %total number of steps
time_step = 1*10^(-9);

% Determination of initial position and velocity
[rvec10,rvec20,v10,v20] = initial_r_v(w0,R_Condon,mbeta);

% from lab frame to center of mass frame
[Rc0,Dc0,vc0,v120] = LabF2CMSF(rvec10,rvec20,v10,v20);

%% simulation
rvec1 = zeros(N,3);
rvec2 = zeros(N,3);
v1 = zeros(N,3);
v2 = zeros(N,3);
Rc = zeros(N,3);
Dc = zeros(N,3);
vc = zeros(N,3);
v12 = zeros(N,3);
E1 = zeros(N,1);
E2 = zeros(N,1);
time_axis = zeros(N,1);

rvec1(1,:) = rvec10;
rvec2(1,:) = rvec20;
v1(1,:) = v10;
v2(1,:) = v20;
Rc(1,:) = Rc0;
Dc(1,:) = Dc0;
vc(1,:) = vc0;
v12(1,:) = v120;
E1(1) = 1/2*m*dot(v10,v10) + U(rvec10);
E2(1) = 1/2*m*dot(v20,v20) + U(rvec20);

i = 2;
while i<=N
    [v1(i-1,:),time_space1] = Cooling(gamma_p,v1(i-1,:),vp);
    [v2(i-1,:),time_space2] = Cooling(gamma_p,v2(i-1,:),vp);
    time_space = max([time_space1,time_space2]); 
    N_cool = floor(time_space/time_step);
    i0 = i;
    while i <= (i0+N_cool)
        [rvec1(i,:),v1(i,:)] = Velocity_Verlet(F,m,time_step,rvec1(i-1,:),v1(i-1,:),kb,TD);
        [rvec2(i,:),v2(i,:)] = Velocity_Verlet(F,m,time_step,rvec2(i-1,:),v2(i-1,:),kb,TD);
        [Rc(i,:),Dc(i,:),vc(i,:),v12(i,:)] = LabF2CMSF(rvec1(i,:),rvec2(i,:),v1(i,:),v2(i,:));
        time_axis(i) = time_axis(i-1) + time_step;
        E1(i) = Energy(rvec1(i,:),v1(i,:),m,U,kb,TD);
        E2(i) = Energy(rvec2(i,:),v2(i,:),m,U,kb,TD);
        if abs(norm(rvec1(i,:) - rvec2(i,:)) - R_Condon)<=Reaction_D
            %transition judgement
            vr = Relative_radial_speed(rvec1(i,:),rvec2(i,:),v1(i,:),v2(i,:));
            Pe = 1 - P_LZ(vr);
            Random_excitation = rand();
            if Random_excitation <= Pe
                %S+P calculation
                for j = i:N-1
                    P_survival = exp(-2*Gamma*time_step*(j-i)); %survival probability
                    Random_survival = rand();
                    if Random_survival > P_survival
                        break
                    end
                    [Dc(j+1,:),v12(j+1,:)] = Velocity_Verlet(Fmol,mu,time_step,Dc(j,:),v12(j,:),kb,TD);
                    Rc(j+1,:) = Rc(j,:) + vc(j,:)*time_step;
                    vc(j+1,:) = vc(j,:);                                    %motion of center of mass
                    [rvec1(j+1,:),rvec2(j+1,:),v1(j+1,:),v2(j+1,:)] = CMSF2LabF(Rc(j+1,:),Dc(j+1,:),vc(j+1,:),v12(j+1,:));
                    time_axis(j+1) = time_axis(j) + time_step;
                end
                i = j;
                E1(i) = Energy(rvec1(i,:),v1(i,:),m,U,kb,TD);
                E2(i) = Energy(rvec2(i,:),v2(i,:),m,U,kb,TD);
                Final_atom = Loss_judgement(E1(i),E2(i));
            end
            if Final_atom == 0 || Final_atom == 1
                break
            end
        end
        i = i + 1;
    end
    i = i + 1;
    if Final_atom == 0 || Final_atom == 1
        break
    end
end
Final_atom = Loss_judgement(E1(i-1),E2(i-1));
