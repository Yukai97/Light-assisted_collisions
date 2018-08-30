function [E] = Energy(rvec,v,m,U,kb,TD)

E = 1/2*m*dot(v,v)/(kb*TD) + U(rvec);


end

