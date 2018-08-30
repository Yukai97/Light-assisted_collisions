function [rvec,v,time_axis,Final_atom] = OnebodyLCal(P_onebodyL,F,U,m,time_step,rvec,v,i0,N_t,time_axis,kb,TD)

EK = Kinetic_Energy(v(i0,:),m,kb,TD);
E = EK + U(rvec(i0,:));
if E>0
    Final_atom = 0;
    return 
end
for i = i0:N_t-1
    [rvec(i+1,:),v(i+1,:)] = Velocity_Verlet(F,m,time_step,rvec(i,:),v(i,:),kb,TD);
    time_axis(i+1) = time_axis(i) + time_step;
    P_oneL = P_onebodyL(time_axis(i+1));
    Random_oneL = rand(1);
    if Random_oneL > P_oneL
        Final_atom = 0;
        break;
    end
end
if i>=N_t
    Final_atom = 1;
end

end

