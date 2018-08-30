function [EK] = Kinetic_Energy(v,m,kb,TD)

EK = 1/2*m*dot(v,v)/(kb*TD);


end

