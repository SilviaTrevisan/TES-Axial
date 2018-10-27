function [c,f,s] = pde_standstill_cycles(x,t,u,DuDx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global Ks_eff rho_s c_s void rho_a c_a D A T_inf
global U_stand
c = [1; 1];
f = [0; Ks_eff/(rho_s*c_s*(1-void))] .* DuDx;
s = [(((U_stand*D*pi)/(rho_a*c_a*void*A))*(T_inf-u(1))); 0];

end

