clear;
clc;

%% parameters
h = 6.626*10^(-34);       %Planck constant, unit: J*s
hbar = h/(2*pi);
g = 9.8;    %gravitational accelaration, unit:m/s^2
kb = 1.38*10^(-23);      %Boltzmann constant, unit: J/K

m = 173.9388*1.6605*10^(-27);   %Yb-174, unit: kg
mu = m/2;             %reduced mass
Gamma = 2*pi*29.1*10^6;          %natural linewidth, unit:Hz
k = 2*pi/(398.9*10^(-9));           %unit: m^(-1)
TD = hbar*Gamma/(2*kb);           %unit: K
gamma_onebodyL = 1;               %one-body loss rate
P_onebodyL = @(t)exp(-gamma_onebodyL*t);

T0 = 10*10^(-6);     %initial temperature, unit: K
mbeta = m/(kb * T0);        %m * beta in Maxwell distribution, unit: kg/J

Delta_Collision = 50*10^6;       %detuning of collision light, unit: Hz
s_collision = 0.05;
Rabi_collision_square = s_collision * Gamma^2/2;    %The square of Rabi frequency

U0 = hbar*70*10^6/(kb*TD);         %trap depth, unit:k_b * T_D
w0 = 0.6 * 10^(-6);   %waist of trap, unit: m
lambda_trap = 532 * 10^(-9); %wavelength of trap laser, unit:m
z_R = pi*w0^2/lambda_trap;      %Rayleigh length, unit:m
U =@(rvec) -U0/(1 + (rvec(3)/z_R)^2) * exp(-2*(rvec(1)^2 + rvec(2)^2)/(w0^2*(1 + rvec(3)^2/z_R^2)));     % trap function
Fx =@(rvec) -4*U0*exp(-2*(rvec(1)^2 + rvec(2)^2)/(w0^2*(1 + rvec(3)^2/z_R^2)))/(w0^2*(1 + (rvec(3)/z_R)^2))*rvec(1);     
Fy =@(rvec) -4*U0*exp(-2*(rvec(1)^2 + rvec(2)^2)/(w0^2*(1 + rvec(3)^2/z_R^2)))/(w0^2*(1 + (rvec(3)/z_R)^2))*rvec(2);
Fz =@(rvec) Trap_Fz(rvec,w0,U0,z_R);    
F = {Fx,Fy,Fz};                    %force corresponding to trap

C3 = 3/2*hbar/(kb*TD)*Gamma/k^3;   %molecular potential constant, unit: kb*TD/m^3
Umol =@(rvec) C3/norm(rvec)^3;
Fmolx = @(rvec) 3*C3/norm(rvec)^5*rvec(1);
Fmoly = @(rvec) 3*C3/norm(rvec)^5*rvec(2);
Fmolz = @(rvec) 3*C3/norm(rvec)^5*rvec(3);       
Fmol = {Fmolx,Fmoly,Fmolz};               %force corresponding to molecular potentials
R_Condon = (C3/(hbar*Delta_Collision)*kb*TD)^(1/3);     %Condon point
Reaction_D = 1/100 * R_Condon;
P_LZ = @(vr) exp(-2*pi*(hbar*Rabi_collision_square/(kb*TD))*(R_Condon^4/C3)/(3*vr));

N_t = 100000;              %total number of steps
time_step = 1*10^(-7);
N_sim = 500;           %number of simulations
Final_atom = 2*ones(1,N_sim);
Loss_time = zeros(1,N_sim);

% Determination of initial position and velocity
[rvec10,rvec20,v10,v20] = initial_r_v(w0,R_Condon,mbeta);

% from lab frame to center of mass frame
[Rc0,Dc0,vc0,v120] = LabF2CMSF(rvec10,rvec20,v10,v20);

%% simulation
tic
for l = 1:N_sim
    rvec1 = zeros(N_t,3);
    rvec2 = zeros(N_t,3);
    v1 = zeros(N_t,3);
    v2 = zeros(N_t,3);
    Rc = zeros(N_t,3);
    Dc = zeros(N_t,3);
    vc = zeros(N_t,3);
    v12 = zeros(N_t,3);
    EK1 = zeros(N_t,1);
    EK2 = zeros(N_t,1);
    time_axis = zeros(N_t,1);

    rvec1(1,:) = rvec10;
    rvec2(1,:) = rvec20;
    v1(1,:) = v10;
    v2(1,:) = v20;
    Rc(1,:) = Rc0;
    Dc(1,:) = Dc0;
    vc(1,:) = vc0;
    v12(1,:) = v120;
    EK1(1) = Kinetic_Energy(v10,m,kb,TD);
    EK2(1) = Kinetic_Energy(v20,m,kb,TD);

    i = 2;
    while i<=N_t
        [rvec1(i,:),v1(i,:)] = Velocity_Verlet(F,m,time_step,rvec1(i-1,:),v1(i-1,:),kb,TD);
        [rvec2(i,:),v2(i,:)] = Velocity_Verlet(F,m,time_step,rvec2(i-1,:),v2(i-1,:),kb,TD);
        [Rc(i,:),Dc(i,:),vc(i,:),v12(i,:)] = LabF2CMSF(rvec1(i,:),rvec2(i,:),v1(i,:),v2(i,:));
        time_axis(i) = time_axis(i-1) + time_step;
        EK1(i) = Kinetic_Energy(v1(i,:),m,kb,TD);
        EK2(i) = Kinetic_Energy(v2(i,:),m,kb,TD);
        P_oneL1 = P_onebodyL(time_axis(i));
        P_oneL2 = P_onebodyL(time_axis(i));
        Random_oneL1 = rand(1);
        Random_oneL2 = rand(1); 
        if Random_oneL1 > P_oneL1 && Random_oneL2 <= P_oneL1
            Final_atom = 1;
            
            break;
        elseif Random_oneL1 <= P_oneL1 && Random_oneL2 > P_oneL2
            Final_atom = 1;
            continue;
        elseif Random_oneL1 > P_oneL1 && Random_oneL2 > P_oneL2
            Final_atom = 0;
            break;
        else
            if abs(norm(rvec1(i,:) - rvec2(i,:)) - R_Condon)<=Reaction_D
                %transition judgement
                vr = Relative_radial_speed(rvec1(i,:),rvec2(i,:),v1(i,:),v2(i,:));
                Pe = 1 - P_LZ(vr);
                Random_excitation = rand(1);
                if Random_excitation <= Pe
                    %S+P calculation
                    for j = i:N_t-1
                        P_survival = exp(-2*Gamma*time_step*(j-i)); %survival probability
                        Random_survival = rand(1);
                        if Random_survival > P_survival
                            break
                        end
                        [Dc(j+1,:),v12(j+1,:)] = Velocity_Verlet(Fmol,mu,time_step,Dc(j,:),v12(j,:),kb,TD);
                        Rc(j+1,:) = Rc(j,:) + vc(j,:)*time_step;
                        vc(j+1,:) = vc(j,:);                                    %motion of center of mass
                        [rvec1(j+1,:),rvec2(j+1,:),v1(j+1,:),v2(j+1,:)] = CMSF2LabF(Rc(j+1,:),Dc(j+1,:),vc(j+1,:),v12(j+1,:));
                        EK1(j+1) = Kinetic_Energy(v1(j+1,:),m,kb,TD);
                        EK2(j+1) = Kinetic_Energy(v2(j+1,:),m,kb,TD);
                        time_axis(j+1) = time_axis(j) + time_step;
                    end
                    i = j;
                    Final_atom(l) = RELoss_judgement(EK1(i),EK2(i),rvec1(i,:),rvec2(i,:),U);
                    if Final_atom(l) == 0 || Final_atom(l) == 1                   
                        break
                    end
                end
            end
        end
        i = i + 1;
    end
    if i>N_t
        Final_atom(l) = RELoss_judgement(EK1(i-1),EK2(i-1),rvec1,rvec2,U);
    end
    Loss_time(l) = i;
end
toc
%% data processing
% Distance = sqrt(sum(abs(Dc).^2,2));
% figure(1);
% plot(Distance);
% title('Distance between two atoms');
% 
% vspeed1 = sqrt(sum(abs(v1).^2,2));
% figure(2);
% plot(vspeed1);
% title('Speed of atom 1');
% 
% vspeed2 = sqrt(sum(abs(v2).^2,2));
% figure(3);
% plot(vspeed2);
% title('Speed of atom 2');
% 
% figure(4); 
% comet3(rvec1(:,1),rvec1(:,2),rvec1(:,3));
% axis([-2*w0 2*w0 -2*w0 2*w0 -z_R z_R]);
% hold on;
% comet3(rvec2(:,1),rvec2(:,2),rvec2(:,3));
% hold off;