clear;
clc;

syms y(t);
ode = diff(y,t) + 1 + 0.2*y + 1000*y*(y-1) == 0;
cond = y(0) == 100;
ySol(t) = dsolve(ode,cond);