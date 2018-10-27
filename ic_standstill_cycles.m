function [ u0 ] = ic_standstill_cycles(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
global Tf_finaldisch Ts_finaldisch xmesh

% dischTf_IC=flip(Tf_finalch);
% dischTs_IC=flip(Ts_finalch);

p_Tf=polyfit(xmesh,Tf_finaldisch,4);
p_Ts=polyfit(xmesh,Ts_finaldisch,4);

u0 = [p_Tf(1)*(x^4)+p_Tf(2)*(x^3)+p_Tf(3)*(x^2)+p_Tf(4)*(x)+p_Tf(5); p_Ts(1)*(x^4)+p_Ts(2)*(x^3)+p_Ts(3)*(x^2)+p_Ts(4)*(x)+p_Ts(5)];
    
end

