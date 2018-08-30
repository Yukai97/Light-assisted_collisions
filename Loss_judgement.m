function [Final_atom] = Loss_judgement(E1,E2)

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

