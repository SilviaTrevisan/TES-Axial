function [pl,ql,pr,qr] = bc_discharging_cycles(xl,ul,xr,ur,t)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global T_inf Tdisch

pl = [0; 0];
ql = [1; 1];
pr = [ur(1)-Tdisch; 0];
qr = [0; 1];

end

