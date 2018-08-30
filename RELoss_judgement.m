function [Final_atom] = RELoss_judgement(EK1,EK2,rvec1,rvec2,U)

U1 = U(rvec1);
U2 = U(rvec2);
E1 = EK1 + U1;
E2 = EK2 + U2;
if E1 > 0 && E2 <= 0 
    Final_atom = 1;
elseif E1 <= 0 && E2 > 0 
    Final_atom = 1;
elseif E1 > 0 && E2 > 0
    Final_atom = 0;
else
    Final_atom = 2;
end


end

