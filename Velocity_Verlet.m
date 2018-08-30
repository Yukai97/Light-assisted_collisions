function [rvec_new,v_new] = Velocity_Verlet(F,m,time_step,rvec,v,kb,TD)
%Velocity Verlet algorithm to calculate the positions and velocities of atoms

a = [F{1}(rvec),F{2}(rvec),F{3}(rvec)]/m*kb*TD;

rvec_new = rvec + v*time_step + 1/2*a*time_step^2;

a_new = [F{1}(rvec_new),F{2}(rvec_new),F{3}(rvec_new)]/m*kb*TD;

v_new = v + 1/2*(a_new + a)*time_step;
end

