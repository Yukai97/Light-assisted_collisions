function [rvec1,rvec2,v1,v2] = CMSF2LabF(Rc,Dc,vc,v12)
% From center of mass frame to lab frame
rvec1 = Rc + 1/2*Dc;
rvec2 = Rc - 1/2*Dc;
v1 = vc + 1/2*v12;
v2 = vc - 1/2*v12;

end

