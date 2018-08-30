function [vr] = Relative_radial_speed(rvec1,rvec2,v1,v2)
%get ralative speed in radial direction
D = rvec1 - rvec2;
e_D = D/norm(D);
vr = abs(dot((v1 - v2),e_D));
end

