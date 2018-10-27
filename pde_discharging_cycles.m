function [c,f,s] = pde_discharging_cycles(x,t,u,DuDx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global rho_s c_s void G h_v U_ext D A T_inf dp rock eps_rock
%global Ks_eff rho_a c_a

% Ks_eff = K_eff_s(u(1)); %[W/m*K]
[k_f,k_s] =  Kunii_Smith (u(1), rock, eps_rock, dp, void);
rho_a = Air_Density(u(1));
c_a = Air_Specific_Heat(u(1));

c = [1; 1];
f = [k_f/(rho_a*c_a*(1-void)); k_s/(rho_s*c_s*(1-void))] .* DuDx;
s = [G/(rho_a*void)*DuDx(1)+(h_v/(rho_a*c_a*void))*(u(2)-u(1))+(((U_ext*D*pi)/(rho_a*c_a*void*A))*(T_inf-u(1))); (h_v/(rho_s*c_s*(1-void)))*(u(1)-u(2))];

end

