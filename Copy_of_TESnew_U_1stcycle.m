clear all
close all

global Ks_eff rho_s c_s void rho_a G h_v c_a U_ext D A T_inf
global Tf_finalch Ts_finalch xmesh
global T_in Tdisch
global Tf_finaldisch Ts_finaldisch
global p_Keff_f p_Keff_s
% global U_stand

% DATA: BASED ON B8 - 
% Fig. 3. Fluid temperature distribution across the storage height at the end of the
% charging period for a series of 6-h charge and 6-h discharge cycles. 
% Baseline parameters:
% H = 3 m, d = 0.02 m, G = 0.225 kg/m2 s, Tin = 823 K, solid material: steatite.

T_in=900+273; %[K]=550°C
Tdisch = 500+273; %[K] discharge T (=20°C)
Tavg = (T_in+Tdisch)/2; %[K] average T
T_inf=293; %[K]

G_real=0.058; %[kg/s*m2] mass flow rate per unit area (A cross section)
% when wall channeling is neglected (D/d_p>40) G_real should be used.
%G=G_real-G_real*0.15; %[kg/s*m2] effective mass flow rate per unit area, it allows to take into account wall wetting phenomena. 15% chosen through exp
G=G_real; %[kg/s*m2] no wall channeling as D/d_p=40

H=1.553; %[m] storage height
D=1.553; %[m] storage diameter
A=(pi*D^2)/4; %[m2]=pi*D^2/4 storage cross section area
void=0.4; %void fraction
d_p=0.02; %[m] praticle diameter
void_dp_D = 0.375+0.17*(d_p/D)+0.39*((d_p/D)^2); %void fraction as f(dp/D)
% void=void_dp_D;
rho_s=2940; %[kg/m3]
c_s=1040; %[J/kg*K]
K_s=2.5; %[W/m*K]
U_ext=0.678; %[W/m2*K] heat transfer coeff through external wall of the storage during working


rho_a=1.0472; %[kg/m3]
c_a=1140; %[J/kg*K]
K_a=0.068; %[W/m*K]
vi=4.1561*10^(-5); %[m2/s]
% rho_a=Air_Density(Tavg); %[kg/m3]
% c_a=Air_Specific_Heat(Tavg);  %[J/kg*K]
% K_a=Air_Conductivity(Tavg); %[W/m*K]
% vi=Air_Dynamic_Viscosity(Tavg); % [kg/ms]
Pr=0.68;

alpha_p=(700/(6*(1-void)))*(G^0.76)*(d_p^0.24);
h_v=(6*(1-void)*alpha_p)/(d_p); %[W/K m3] volumetric heat transfer coeff

% Ks_eff=((void*K_a^-1)+((1-void)*K_s^-1))^-1; %[W/m*K] effective solid conductivity (weighted mean between solid and fluid)
Ks_eff=void*K_a+(1-void)*K_s; %[W/m*K] effective solid conductivity (weighted mean between solid and fluid)
% Effective conductivity of the bed fluid and solid: condutivity, radiation (NO MIXING) 
% loaded_keff = 'keff_Kunii_Smith_Basalt.mat';
% load(loaded_keff,'p_lambda_eff_f','p_lambda_eff_s');
% p_Keff_f = p_lambda_eff_f; p_Keff_s = p_lambda_eff_s;

%for Ansys input
permeability=(d_p^2/150)*(void^3/((1-void)^2));
C2=(3.5/d_p)*((1-void)/void^3);  %inertial_loss_coeff

%% PDE solution___ CHARGING___

m=0;
tmax=12*3600; %6h of charging cycle
% dx = 0.01;
% xmesh = 0:dx:H; %[m]
% [~,Nr] = size(xmesh); % number of steps in space
n_xmesh=56; % 1 step = 0.01 m
n_tstep=721; % 1 step = 1 min
xmesh=linspace(0,H,n_xmesh);
tstep=linspace(0,tmax,n_tstep); % 1 step = 1 min

sol = pdepe(m,@pde_charge_cycles,@ic_charge_cycles,@bc_charge_cycles,xmesh,tstep);
Tf = sol(:,:,1); Ts = sol(:,:,2);
Tf_finalch=Tf(end,:); Ts_finalch=Ts(end,:);
T_ic=[Tf_finalch; Ts_finalch];

%% Graphs charging

figure 
plot(xmesh,Tf(11,:))
hold on
plot(xmesh,Tf(61,:))
plot(xmesh,Tf(181,:))
plot(xmesh,Tf(301,:))
plot(xmesh,Tf(end,:))
% plot(xmesh,Ts(20,:),':')
% plot(xmesh,Ts(50,:),':')
% plot(xmesh,Ts(80,:),':')
% plot(xmesh,Ts(end,:),':')
hold off
xlabel('x');
ylabel('Tf');
title('CHARGING');
legend('Tf @t=10 min', 'Tf @t=60 min', 'Tf @t=180 min', 'Tf @t=300 min', 'Tf end')

%% Pressure drop 
%correlation from B8, valid for 1<Re<10^4
Re=(G*d_p)/(vi*rho_a);
Pr_formula=c_a*vi*rho_a*K_a;
%G is the effective mass flow rate per unit area, reduced to take into
%account wall channeling. Therefore, Re and pressure drop are lower. 
coeff_dp=1.75*((1-void)/(void^3))+150*((1-void)/(void^3))*((rho_a*vi)/(G*d_p));
delta_pressure=((H*G^2)/(rho_a*d_p))*coeff_dp; %[Pa] pressure drop along the packed bed

%% PDE solution____DISCHARGING______

% sol_disch = pdepe(m,@TESnewU_pde_cycles,@TESnewU_ic_discharging_cycles,@TESnewU_bc_discharging_cycles,xmesh,tstep);

sol_disch = pdepe(m,@pde_discharging_cycles,@ic_discharging2_cycles,@bc_discharging2_cycles,xmesh,tstep);

Tf_disch = sol_disch(:,:,1); Ts_disch = sol_disch(:,:,2);

Tf_finaldisch=Tf_disch(end,:); Ts_finaldisch=Ts_disch(end,:);

%% Graphs discharging
figure 
plot(xmesh,Tf_disch(11,:))
hold on
plot(xmesh,Tf_disch(61,:))
plot(xmesh,Tf_disch(181,:))
plot(xmesh,Tf_disch(301,:))
plot(xmesh,Tf_disch(end,:))
% plot(xmesh,Ts_disch(20,:),':')
% plot(xmesh,Ts_disch(50,:),':')
% plot(xmesh,Ts_disch(80,:),':')
% plot(xmesh,Ts_disch(end,:),':')
hold off
xlabel('x');
ylabel('Tf');
title('DISCHARGING');
legend('Tf @t=10 min', 'Tf @t=60 min', 'Tf @t=180 min', 'Tf @t=300 min', 'Tf @t=480 min')

%% STORAGE RATE OF CHARGE/DISCHARGE:
dt=(tmax/60)/(n_tstep-1); %[min] step length over time
Rate_Charge=diff(Ts)./dt; %[°C/min]
Rate_Discharge=diff(Ts_disch)./dt; %[°C/min]
time_rate=tstep(1:end-1);

% figure 
% surf(xmesh,time_rate,Rate_Charge);
% title('Rate of charge');
% xlabel('x');
% ylabel('t');
% figure 
% surf(xmesh,time_rate,Rate_Discharge);
% title('Rate of discharge');
% xlabel('x');
% ylabel('t');

%% PERFORMANCE INDICATORS: attention when properties depends on T

% E_input, Energy released (lost) by the fluid during the charge phase:
% deltaT_input=Tf(:,1)-Tf(:,end); deltaT_input=T_in-T_inf; E_input=A*G_real*c_a*deltaT_input*tmax; %[J]
DeltaT_input=T_in-Tf(:,end); %[K] delta T between inlet and outlet of the storage for each timestep
E_in_tstep=A*G*c_a*dt*60.*DeltaT_input; %[J] Einput in each timestep
E_input=sum(E_in_tstep); %[J]
E_input_kWh=2.7778e-7*E_input; %[kWh]

% E_output, Energy recovered by the fluid during the discharging phase:
% deltaT_output=Tf_disch(:,1)-T_inf; E_output=A*G_real*c_a*(trapz(deltaT_output)*(tmax/180)); %[J] ?????
DeltaT_outflow=T_inf-Tf_disch(:,1); %[K] delta T between inlet and outlet of the storage for each timestep
E_out_tstep=A*G*c_a*dt*60.*DeltaT_outflow; %[J] Eoutflow in each timestep
E_outflow=sum(E_out_tstep); %[J]
E_outflow_kWh=2.7778e-7*E_outflow; %[kWh]

% Energy stored in the filler material:
DeltaT_stored=Ts(end,:)-T_inf; %[K] delta T between atorage inlet and outlet at the end of charge
E_stored=rho_s*A*c_s*(1-void)*(trapz(DeltaT_stored)*(H/(n_xmesh-1))); %[J]
E_stored_max=rho_s*A*c_s*(1-void)*(T_in-T_inf)*H;  %[J]
E_stored_kWh=2.7778e-7*E_stored; %[kWh]

% Energy of pumping (isothermal compressor):
E_pump=(A*G*delta_pressure*tmax)/rho_a; %[J]

% Performance indicators:
eta_charging=E_stored/(E_input+E_pump);
eta_discharging=abs(E_outflow)/(E_stored+E_pump);
eta_overall=abs(E_outflow)/(E_input+E_pump+E_pump);
CR=E_stored/E_stored_max; % Capacity ratio

%% PDE solution____STANDSTILL____ Check Standstill solution via B4

% alpha_i=(K_a/d_p)*(2.58*(Re^(1/3))+0.094*(Re^0.8)*(Pr^0.4));
% U_stand=((1/U_ext)-(1/alpha_i))^-1; % [W/m2*K] heat transfer coeff through external wall of the storage durign standstill
% % why U_stand > U_ext? Is it feasible??
% 
% tmax_standstill=12*3600; %12h of standstill
% tstep_standstill=linspace(0,tmax_standstill,360); % 1 step = 2 min
% 
% sol_stand = pdepe(m,@pde_standstill_cycles,@ic_standstill_cycles,@bc_standstill_cycles,xmesh,tstep_standstill);
% 
% Tf_stand = sol_stand(:,:,1);
% Ts_stand = sol_stand(:,:,2);
% 
% figure 
% plot(xmesh,Tf_stand(30,:),'r')
% hold on
% plot(xmesh,Tf_stand(90,:),'g')
% plot(xmesh,Tf_stand(240,:),'b')
% plot(xmesh,Tf_stand(end,:),'c')
% plot(xmesh,Ts_stand(30,:),'r:')
% plot(xmesh,Ts_stand(90,:),'g:')
% plot(xmesh,Ts_stand(240,:),'b:')
% plot(xmesh,Ts_stand(end,:),'c:')
% hold off
% xlabel('x');
% ylabel('Tf and Ts');
% title('STANDSTILL');
% legend('Tf @t=1h', 'Tf @t=3h', 'Tf @t=8h', 'Tf @t=end', 'Ts @t=1h', 'Ts @t=3h', 'Ts @t=8h', 'Ts @t=end')
