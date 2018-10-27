clear all
close all

global Ks_eff rho_s c_s void rho_a G h_v c_a U_ext D A T_inf
global Tf_final Ts_final xmesh
global T_in Tdisch
global p_Keff_f p_Keff_s p
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

m = 0.1101; %[kg/s] air mass flow rate
G_real=m/A; %[kg/s*m2] mass flow rate per unit area (A cross section)
% when wall channeling is neglected (D/d_p>40) G_real should be used.
%G=G_real-G_real*0.15; %[kg/s*m2] effective mass flow rate per unit area, it allows to take into account wall wetting phenomena. 15% chosen through exp
G=G_real; %[kg/s*m2] no wall channeling as D/d_p=40

p = 300000;
rho_a=Air_Density(Tavg); %[kg/m3]
c_a=Air_Specific_Heat(Tavg);  %[J/kg*K]
K_a=Air_Conductivity(Tavg); %[W/m*K]
vi=Air_Dynamic_Viscosity(Tavg); % [kg/ms]
Pr=0.68;

alpha_p=(700/(6*(1-void)))*(G^0.76)*(d_p^0.24); %[W/K m2] heat transfer coeff air/bed
h_v=(6*(1-void)*alpha_p)/(d_p); %[W/K m3] volumetric heat transfer coeff

% Ks_eff=((void*K_a^-1)+((1-void)*K_s^-1))^-1; %[W/m*K] effective solid conductivity (weighted mean between solid and fluid)
Ks_eff=void*K_a+(1-void)*K_s; %[W/m*K] effective solid conductivity (weighted mean between solid and fluid)
% Effective conductivity of the bed fluid and solid: condutivity, radiation (NO MIXING) 
loaded_keff = 'keff_Kunii_Smith_Basalt.mat';
load(loaded_keff,'p_lambda_eff_f','p_lambda_eff_s');
p_Keff_f = p_lambda_eff_f; p_Keff_s = p_lambda_eff_s;

%for Ansys input
permeability=(d_p^2/150)*(void^3/((1-void)^2));
C2=(3.5/d_p)*((1-void)/void^3);  %inertial_loss_coeff

%% PDE solution___ CHARGING___

m=0;
tmax=12*3600; %6h of charging cycle
% dx = 0.1;
% xmesh = 0:dx:H; %[m]
% [~,Nr] = size(xmesh); % number of steps in space
n_xmesh=78; % 1 step = 0.01 m
n_tstep=721; % 1 step = 1 min
dt = tmax/(n_tstep-1);
dx = H/(n_xmesh-1);
xmesh=linspace(0,H,n_xmesh);
tstep=linspace(0,tmax,n_tstep); % 1 step = 1 min
xmesh_adim = xmesh./H;

sol = pdepe(m,@pde_charge_cycles,@ic_charge_cycles,@bc_charge_cycles,xmesh,tstep);
Tf_ch = sol(:,:,1); Ts_ch = sol(:,:,2);
Tf_final=Tf_ch(end,:); Ts_final=Ts_ch(end,:);

%% Graphs charging

figure(1)
plot(xmesh,Tf_ch(11,:), xmesh,Tf_ch(61,:), xmesh,Tf_ch(181,:), xmesh,Tf_ch(301,:), xmesh,Tf_ch(end,:) )
xlabel('x'), ylabel('Tf'), title('CHARGING'),
legend('Tf @t=10 min', 'Tf @t=60 min', 'Tf @t=180 min', 'Tf @t=300 min', 'Tf end')

%% Pressure drop 
%correlation from B8, valid for 1<Re<10^4
Re=(G*d_p)/vi;
v_sup = G/rho_a;
v_int = v_sup/void;
Pr_formula=c_a*vi/K_a;
coeff_dp=1.75*((1-void)/(void^3))+150*(((1-void)^2)/(void^3))*((vi)/(G*d_p));
delta_pressure=((H*G^2)/(rho_a*d_p))*coeff_dp; %[Pa] pressure drop along the packed bed

%% PDE solution____DISCHARGING______

sol_disch = pdepe(m,@pde_discharging_cycles,@ic_discharging2_cycles,@bc_discharging2_cycles,xmesh,tstep);
Tf_disch = sol_disch(:,:,1); Ts_disch = sol_disch(:,:,2);
Tf_final=Tf_disch(end,:); Ts_final=Ts_disch(end,:);

%% Graphs discharging
figure 
plot(xmesh,Tf_disch(11,:), xmesh,Tf_disch(61,:) ,xmesh,Tf_disch(181,:), xmesh,Tf_disch(301,:), xmesh,Tf_disch(end,:))
xlabel('x'), ylabel('Tf'), title('DISCHARGING'),
legend('Tf @t=10 min', 'Tf @t=60 min', 'Tf @t=180 min', 'Tf @t=300 min', 'Tf @t=480 min')


%% PDE solution _____CYCLES____

for n=1:9
    sol = pdepe(m,@pde_charge_cycles,@ic_discharging2_cycles,@bc_charge_cycles,xmesh,tstep);
    Tf_ch_cycles(:,:,n) = sol(:,:,1); Ts_ch_cycles(:,:,n) = sol(:,:,2);
    Tf_final=Tf_ch_cycles(end,:,n); Ts_final=Ts_ch_cycles(end,:,n);

    sol_disch = pdepe(m,@pde_discharging_cycles,@ic_discharging2_cycles,@bc_discharging2_cycles,xmesh,tstep);
    Tf_disch_cycles(:,:,n) = sol_disch(:,:,1); Ts_disch_cycles(:,:,n) = sol_disch(:,:,2);
    Tf_final=Tf_disch_cycles(end,:,n); Ts_final=Ts_disch_cycles(end,:,n);
end
%% Graphs cycles
figure(3)
plot(xmesh,Tf_ch(361,:),'b--', xmesh,Tf_ch(end,:),'b', xmesh,Tf_ch_cycles(361,:,2),'r--', xmesh,Tf_ch_cycles(end,:,2),'r', xmesh,Tf_ch_cycles(361,:,end),'k--', xmesh,Tf_ch_cycles(end,:,end),'k')
xlabel('x'), ylabel('Tf'), title('Charge - cycles'),
legend('1st cycle: Tf @t=360 min', 'Tf @t=720 min','3rd cycle: Tf @t=360 min', 'Tf @t=720 min','10th cycle: Tf @t=360 min', 'Tf @t=720 min')

%% External Radiation
eps = 0.1;
sigma = 5.67*10^-8;
DT_wall = 300;
Qrad_lat = eps*sigma*pi*D*trapz(xmesh,((Tf_ch-DT_wall).^4-T_inf^4)');
Erad_lat = sum(Qrad_lat.*dt);
Qrad_top = eps*sigma*pi*(D^2/4)*((Tf_ch(:,1)-DT_wall).^4-T_inf^4);
Erad_top = sum(Qrad_top.*dt);

Erad_tot = Erad_lat+Erad_top;
%% PERFORMANCE INDICATORS: attention when properties depends on T

% E_input, Energy released (lost) by the fluid during the charge phase:
E_in_tstep=A*G*c_a*dt*(T_in-Tf_ch(:,end)); %[J] Einput in each timestep
E_input=sum(E_in_tstep); %[J] 
E_input_kWh=2.7778e-7*E_input; %[kWh]

% E_output, Energy recovered by the fluid during the discharging phase:
E_out_tstep=A*G*c_a*dt.*(Tf_disch(:,1)-Tf_disch(:,end)); %[J] Eoutflow in each timestep
E_outflow=sum(E_out_tstep); %[J]
E_outflow_kWh=2.7778e-7*E_outflow; %[kWh]

% Energy stored in the filler material:
E_stored=rho_s*A*c_s*(1-void)*trapz(xmesh, Ts_ch(end,:)-Ts_ch(1,:)); %[J]
E_stored_max=rho_s*A*c_s*(1-void)*(T_in-Tdisch)*H;  %[J]
E_stored_kWh=2.7778e-7*E_stored; %[kWh]

% Energy of pumping (isothermal compressor):
E_pump=(A*G*delta_pressure*tmax)/rho_a; %[J]

% Performance indicators:
eta_charging=E_stored/(E_input+E_pump); eta_discharging=abs(E_outflow)/(E_stored+E_pump);
eta_overall=abs(E_outflow)/(E_input+E_pump+E_pump); CR=E_stored/E_stored_max; % Capacity ratio

cycles = 1:10;
for n=1:9
    % Energy:
    E_in_tstep_cycles(:,n)=A*G*c_a*dt*(T_in-Tf_ch_cycles(:,end,n)); %[J] Einput in each timestep
    E_input_cycles(n)=sum(E_in_tstep_cycles(:,n)); %[J] 
    E_out_tstep_cycles(:,n)=A*G*c_a*dt.*(Tf_disch_cycles(:,1,n)-Tf_disch_cycles(:,end,n)); %[J] Eoutflow in each timestep
    E_outflow_cycles(n)=sum(E_out_tstep_cycles(:,n)); %[J]
    E_stored_cycles(n)=rho_s*A*c_s*(1-void)*trapz(xmesh, Ts_ch_cycles(end,:,n)-Ts_ch_cycles(1,:,n)); %[J]
    % Performance indicators:
    eta_charging_cycles(n)=E_stored_cycles(n)/(E_input_cycles(n)+E_pump); 
    eta_discharging_cycles(n)=abs(E_outflow_cycles(n))/(E_stored_cycles(n)+E_pump);
    eta_overall_cycles(n)=abs(E_outflow_cycles(n))/(E_input_cycles(n)+E_pump+E_pump); 
    CR_cycles(n)=E_stored_cycles(n)/E_stored_max; % Capacity ratio
end

eta_ch = [eta_charging, eta_charging_cycles];
eta_disch = [eta_discharging, eta_discharging_cycles];
eta_ovr = [eta_overall, eta_overall_cycles];
CR_tot = [CR, CR_cycles];

plot(cycles, eta_ch, cycles, eta_disch, cycles, eta_ovr, cycles, CR_tot)
xlabel('cycles'), ylabel('eta [-]'), title('efficiency'), grid on,
legend('charge','discharge','overall','Capacity Ratio')