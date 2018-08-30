function [Fz] = Trap_Fz(rvec,w0,U0,z_R)
%Force on atom of trap in z direction

wz = sqrt(1 + (rvec(3)/z_R)^2)*w0;
Fz1 = 4*U0*rvec(3)*exp(-2*(rvec(1)^2 + rvec(2)^2)/wz^2)*(rvec(1)^2 + rvec(2)^2)/(w0^2*(1 + rvec(3)^2/z_R^2)^3*z_R^2);
Fz2 = -2*U0*rvec(3)*exp(-2*(rvec(1)^2 + rvec(2)^2))/((1+rvec(3)^2/z_R^2)^2*z_R^2);
Fz = Fz1 + Fz2;

end

