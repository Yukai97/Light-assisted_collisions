function [vcooling,time_space] = Cooling(gamma_p,v,vp)

vfc = [v,-v];
vcooling = v;
rate = gamma_p(vfc);
prob_norm = rate/sum(rate);
random_num = rand();
if  0 <= random_num && random_num < prob_norm(1)
    vcooling(1) = vcooling(1) - vp;
    time_space = 1/rate(1);
elseif prob_norm(1) <= random_num && random_num < (prob_norm(1) + prob_norm(2))
    vcooling(1) = vcooling(1) + vp;
    time_space = 1/rate(2);
elseif sum(prob_norm(1:2)) <= random_num && random_num < sum(prob_norm(1:3))
    vcooling(2) = vcooling(2) - vp;
    time_space = 1/rate(3);
elseif sum(prob_norm(1:3)) <= random_num && random_num < sum(prob_norm(1:4))
    vcooling(2) = vcooling(2) + vp;
    time_space = 1/rate(4);
elseif sum(prob_norm(1:4)) <= random_num && random_num < sum(prob_norm(1:5))
    vcooling(3) = vcooling(3) - vp;
    time_space = 1/rate(5);
else 
    vcooling(3) = vcooling(3) + vp;
    time_space = 1/rate(6);
end

end

