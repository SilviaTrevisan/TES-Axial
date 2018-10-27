function [pl,ql,pr,qr] = bc_charge_cycles(xl,ul,xr,ur,t)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global T_in
pl = [ul(1)-T_in; 0];
ql = [0; 1];
pr = [0; 0];
qr = [1; 1];

end

