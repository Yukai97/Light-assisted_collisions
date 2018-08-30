function [Rc,Dc,vc,v12] = LabF2CMSF(rvec1,rvec2,v1,v2)
%From lab frame to center of mass frame

Rc = 1/2*(rvec1 + rvec2);
Dc = rvec1 - rvec2;
vc = 1/2*(v1 + v2);
v12 = v1 - v2;

end

