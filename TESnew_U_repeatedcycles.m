clear all
close all

global Ks_eff rho_s c_s void rho_a G h_v c_a U_ext D A T_inf
global Tf_finalch Ts_finalch xmesh
global T_in
global Tf_finaldisch Ts_finaldisch

% DATA: BASED ON B8 - 
% Fig. 3. Fluid temperature distribution across the storage height at the end of the
% charging period for a series of 6-h charge and 6-h discharge cycles. 
% Baseline parameters:
% H = 3 m, d = 0.02 m, G = 0.225 kg/m2 s, Tin = 823 K, solid material: steatite.

T_in=823; %[K]=550°C
T_inf=293; %[K]

G_real=0.225; %[kg/s*m2] mass flow rate per unit area (A cross section)
% G_real=0.078; %[kg/s*m2] mass flow rate per unit area (A cross section), cosidering the same mass flow rate as in table 1 (B8)
% when wall channeling is neglected (D/d_p>40) G_real should be used.
% G=G_real-G_real*0.15; %[kg/s*m2] effective mass flow rate per unit area, it allows to take into account wall wetting phenomena. 15% chosen through exp
G=G_real; %[kg/s*m2] no wall channeling as D/d_p=40

H=3; %[m] storage height
D=0.8; %[m] storage diameter
A=(pi*D^2)/4; %[m2]=pi*D^2/4 storage cross section area
void=0.4; %void fraction
d_p=0.02; %[m] praticle diameter
rho_s=2680; %[kg/m3]
c_s=1068; %[J/kg*K]
K_s=2.5; %[W/m*K]
U_ext=0.678; %[W/m2*K] heat transfer coeff through external wall of the storage

rho_a=0.6337; %[kg/m3]
c_a=1043.1; %[J/kg*K]
K_a=0.04441; %[W/m*K]
vi=45.846*10^(-6); %[m2/s]
Pr=0.68;

alpha_p=(700/(6*(1-void)))*(G^0.76)*(d_p^0.24);
h_v=(6*(1-void)*alpha_p)/(d_p); %[W/K m3] volumetric heat transfer coeff

% Ks_eff=((void*K_a^-1)+((1-void)*K_s^-1))^-1; %[W/m*K] effective solid conductivity (weighted mean between solid and fluid)
Ks_eff=void*K_a+(1-void)*K_s; %[W/m*K] effective solid conductivity (weighted mean between solid and fluid)

%for Ansys input
permeability=(d_p^2/150)*(void^3/((1-void)^2));
C2=(3.5/d_p)*((1-void)/void^3);  %inertial_loss_coeff

%% PDE solution___ 1st CHARGING___

m=0;
tmax=6*3600; % 6h of charging cycle
n_xmesh=50; % 1 step = 0.06 m
n_tstep=360; % 1 step = 1 min
xmesh=linspace(0,H,n_xmesh);
tstep=linspace(0,tmax,n_tstep); 

sol = pdepe(m,@pde_charge_cycles,@ic_charge_cycles,@bc_charge_cycles,xmesh,tstep);
Tf = sol(:,:,1);
Ts = sol(:,:,2);

Tf_finalch=Tf(end,:);
Ts_finalch=Ts(end,:);
%T_ic (2 x xmesh) vector with the initial condition for the discharging
%cycle (equal to the final condition of the charging cycle)
T_ic=[Tf_finalch; Ts_finalch];

%% Graphs charging

% figure 
% surf(xmesh,tstep,Tf);
% title('Tf');
% xlabel('x');
% ylabel('t');
% 
% figure 
% surf(xmesh,tstep,Ts);
% title('Ts');
% xlabel('x');
% ylabel('t');

figure % T fluid and T solid along the storage lenght with time as parameter
plot(xmesh,Tf(20,:))
hold on
plot(xmesh,Tf(50,:))
plot(xmesh,Tf(80,:))
plot(xmesh,Tf(end,:))
plot(xmesh,Ts(20,:),':')
plot(xmesh,Ts(50,:),':')
plot(xmesh,Ts(80,:),':')
plot(xmesh,Ts(end,:),':')
hold off
xlabel('x');
ylabel('Tf and Ts');
title('CHARGING');
legend('Tf @t=1200s', 'Tf @t=3000s', 'Tf @t=4800s', 'Tf @end charging', 'Ts @t=1200s', 'Ts @t=3000s', 'Ts @t=4800s', 'Ts @end charging')

%% Pressure drop 
%correlation from B8, valid for 1<Re<10^4
Re=(G*d_p)/(vi*rho_a);
%G is the effective mass flow rate per unit area, reduced to take into
%account wall channeling. Therefore, Re and pressure drop are lower. 
coeff_dp=1.75*((1-void)/(void^3))+150*((1-void)/(void^3))*((rho_a*vi)/(G*d_p));
delta_pressure=((H*G^2)/(rho_a*d_p))*coeff_dp; %[Pa] pressure drop along the packed bed

%% PDE solution____CONSECUTIVE DISCHARGE/CHARGE

n_max=20; %number of cycles considered
cycles=1:n_max;
Tf_disch_cycle = zeros(n_tstep,n_xmesh,n_max);
Ts_disch_cycle = zeros(n_tstep,n_xmesh,n_max);
Tf_ch_cycle = zeros(n_tstep,n_xmesh,n_max);
Ts_ch_cycle = zeros(n_tstep,n_xmesh,n_max);

for n = 1:n_max
    
    % DISCHARGE
    % sol_disch = pdepe(m,@TESnewU_pde_cycles,@TESnewU_ic_discharging_cycles,@TESnewU_bc_discharging_cycles,xmesh,tstep);
    sol_disch = pdepe(m,@pde_discharging_cycles,@ic_discharging2_cycles,@bc_discharging2_cycles,xmesh,tstep);

    Tf_disch_cycle(:,:,n) = sol_disch(:,:,1);
    Ts_disch_cycle(:,:,n) = sol_disch(:,:,2);

    Tf_finaldisch=Tf_disch_cycle(end,:,n);
    Ts_finaldisch=Ts_disch_cycle(end,:,n);
    
    % CHARGE
    sol = pdepe(m,@pde_charge_cycles,@ic_charge2_cycles,@bc_charge_cycles,xmesh,tstep);
    
    Tf_ch_cycle(:,:,n) = sol(:,:,1);
    Ts_ch_cycle(:,:,n) = sol(:,:,2);

    Tf_finalch=Tf_ch_cycle(end,:,n);
    Ts_finalch=Ts_ch_cycle(end,:,n);
end
    
%% Graphs cycles
% Graphs from literature paper B8 
cycle1=[0.08370 821.87134;
0.18605 820.01565;
0.28840 817.22413;
0.39761 812.80021;
0.46705 806.60159;
0.53630 797.89973;
0.60220 783.82075;
0.67047 763.52344;
0.73194 740.37482;
0.78481 716.24510;
0.83268 692.20366;
0.87756 668.88704;
0.91944 644.61564;
0.96131 619.73346;
1.00019 596.30139;
1.03907 572.41124;
1.08095 546.91831;
1.12282 522.64690;
1.16470 498.52818;
1.20659 476.24176;
1.25147 453.11601;
1.30234 430.14388;
1.35621 406.99408;
1.41489 385.27015;
1.48196 364.58398;
1.54545 349.16372;
1.61374 334.48268;
1.67891 323.85303;
1.74386 316.45763;
1.80925 309.71123;
1.87843 304.82317;
1.94166 301.83010;
2.02730 298.78963;
2.11422 296.66550;
2.26710 294.55167;
2.33305 294.56195;
2.50175 294.58827;
2.69477 294.05545;
2.79769 294.02875;
2.89362 293.77769;
2.98355 294.16937];

cycle3=[0.08632, 821.87175;
0.25691, 818.15504;
0.39761, 812.80021;
0.46707, 811.28074;
0.53633, 804.91845;
0.64120, 791.64837;
0.72408, 776.71416;
0.80294, 756.06222;
0.93630, 712.51235;
1.04244, 670.59237;
1.12924, 634.26534;
1.20107, 603.28018;
1.26691, 575.50061;
1.32678, 550.92662;
1.46148, 498.57447;
1.57824, 460.11444;
1.64412, 439.81675;
1.68903, 426.99768;
1.76690, 407.77070;
1.84178, 390.68094;
2.01953, 358.64346;
2.11261, 345.40436;
2.21528, 333.73439;
2.30518, 325.19769;
2.39909, 318.08674;
2.50698, 311.69052;
2.59690, 307.42919;
2.69481, 303.88166;
2.79772, 301.04747;
2.89364, 298.92475;
2.97459, 301.07506];

figure 
plot(xmesh,Ts(end,:),'b')
hold on
plot(xmesh,Ts_disch_cycle(end,:,1),'r')
plot(xmesh,Ts_ch_cycle(end,:,1),'g')
plot(cycle1(:,1),cycle1(:,2),'b:')
plot(cycle3(:,1),cycle3(:,2),'g:')
hold off
xlabel('x');
ylabel('Ts [K]');
title('First cycles');
legend('end 1st charge','end 1st discharge','end 2nd charge', 'end 1st charge-paper','end 2nd charge-parer')

figure 
plot(xmesh,Ts(end,:),'b')
hold on
plot(xmesh,Ts_ch_cycle(end,:,1),'r')
plot(xmesh,Ts_ch_cycle(end,:,3),'g')
plot(xmesh,Ts_ch_cycle(end,:,10),'k')
plot(xmesh,Ts_ch_cycle(end,:,20),'c')
hold off
xlabel('x');
ylabel('Ts [K]');
title('Cycles charge - solid T along x at the end of charge phase');
legend('n=0','n=1','n=3','n=10','n=20')

figure 
plot(xmesh,Ts_disch_cycle(end,:,1),'r')
hold on
plot(xmesh,Ts_disch_cycle(end,:,3),'g')
plot(xmesh,Ts_disch_cycle(end,:,10),'k')
plot(xmesh,Ts_disch_cycle(end,:,20),'c')
hold off
xlabel('x');
ylabel('Ts [K]');
title('Cycles discharge - solid T along x at the end of discharge phase');
legend('n=1','n=3','n=10','n=20')

figure 
plot(tstep,Tf_disch_cycle(:,1,1),'r')
hold on
plot(tstep,Tf_disch_cycle(:,1,3),'g')
plot(tstep,Tf_disch_cycle(:,1,10),'k')
plot(tstep,Tf_disch_cycle(:,1,20),'c')
hold off
xlabel('time [s]');
ylabel('Tf [K]');
title('Cycles discharge - solid T in time at storage top during discharge');
legend('n=1','n=3','n=10','n=20')

% STORAGE RATE OF CHARGE/DISCHARGE:
dt=(tmax/60)/(n_tstep-1); %[min] step length over time

% Rate of charge/discharge in cyclic conditions!!
%% PERFORMANCE INDICATORS:

% 1st charging cycle:
% E_input, Energy released (lost) by the fluid during the charge phase:
DeltaT_input_1st=T_in-Tf(:,end); %[K] delta T between inlet and outlet of the storage for each timestep
E_in_tstep_1st=A*G*c_a*dt*60.*DeltaT_input_1st; %[J] Einput in each timestep
E_input_1st=sum(E_in_tstep_1st); %[J]
% Energy stored in the filler material:
DeltaT_stored_1st=Ts(end,:)-T_inf; %[K] delta T between atorage inlet and outlet at the end of charge
E_stored_1st=rho_s*A*c_s*(1-void)*(trapz(DeltaT_stored_1st)*(H/(n_xmesh-1))); %[J]
E_stored_max=rho_s*A*c_s*(1-void)*(T_in-T_inf)*H;  %[J]
% Energy of pumping (isothermal compressor):
E_pump=(A*G*delta_pressure*tmax)/rho_a; %[J]

for n = 1:n_max
    
    % E_output (from 1st discharge cycle), Energy recovered by the fluid during the discharging phase:
    DeltaT_outflow_cycle(:,1,n)=T_inf-Tf_disch_cycle(:,1,n); %[K] delta T between inlet and outlet of the storage for each timestep
    E_out_tstep_cycle(:,1,n)=A*G*c_a*dt*60.*DeltaT_outflow_cycle(:,1,n); %[J] Eoutflow in each timestep
    E_outflow_cycle(n,1)=sum(E_out_tstep_cycle(:,1,n)); %[J]
       
    % E_input (form 2nd charge cycle), Energy released (lost) by the fluid during the charge phase:
    DeltaT_input_cycle(:,1,n)=T_in-Tf_ch_cycle(:,end,n); %[K] delta T between inlet and outlet of the storage for each timestep
    E_in_tstep_cycle(:,1,n)=A*G*c_a*dt*60.*DeltaT_input_cycle(:,1,n); %[J] Einput in each timestep
    E_input_cycle(n,1)=sum(E_in_tstep_cycle(:,1,n)); %[J] column vector each element is a cycle

    % Energy stored in the filler material (form 2nd charge cycle):
    DeltaT_sto_cycle_afterch(1,:,n)=Ts_ch_cycle(end,:,n)-T_inf; %[K] delta T between atorage inlet and outlet at the end of charge
    E_stored_cycle_afterch(n,1)=rho_s*A*c_s*(1-void)*(trapz(DeltaT_sto_cycle_afterch(1,:,n))*(H/(n_xmesh-1))); %[J]

    DeltaT_sto_cycle_beforech(1,:,n)=Ts_disch_cycle(end,:,n)-T_inf; %[K] delta T between atorage inlet and outlet at the end of charge
    E_stored_cycle_beforech(n,1)=rho_s*A*c_s*(1-void)*(trapz(DeltaT_sto_cycle_beforech(1,:,n))*(H/(n_xmesh-1))); %[J]
    
    E_stored_cycle(n,1)=E_stored_cycle_afterch(n,1)-E_stored_cycle_beforech(n,1);
    
end

%Energy in kWh:
E_outflow_cycle_kWh=2.7778e-7*E_outflow_cycle;
E_input_cycle_kWh=2.7778e-7*E_input_cycle;
E_stored_cycle_kWh=2.7778e-7*E_stored_cycle;

% Performance indicators: include storage efficiency if we consider standstill!!
% 1st cycle (1st charge + 1st discharge)
eta_charging_cycle(1,1)=E_stored_1st/(E_input_1st+E_pump);
eta_discharging_cycle(1,1)=abs(E_outflow_cycle(1,1))/(E_stored_1st+E_pump);
eta_overall_cycle(1,1)=abs(E_outflow_cycle(1,1))/(E_input_1st+E_pump+E_pump);
CR_cycle(1,1)=E_stored_1st/E_stored_max; % Capacity ratio

for n = 2:n_max
    % from 2nd cycle (2nd charge + 2nd discharge)
    eta_charging_cycle(n,1)=E_stored_cycle(n-1,1)/(E_input_cycle(n-1,1)+E_pump);
    eta_discharging_cycle(n,1)=abs(E_outflow_cycle(n,1))/(E_stored_cycle(n-1,1)+E_pump);
    eta_overall_cycle(n,1)=abs(E_outflow_cycle(n,1))/(E_input_cycle(n-1,1)+E_pump+E_pump);
    CR_cycle(n,1)=E_stored_cycle(n-1,1)/E_stored_max; % Capacity ratio

end

figure 
plot(cycles,eta_charging_cycle)
hold on
plot(cycles,eta_discharging_cycle)
plot(cycles,eta_overall_cycle)
plot(cycles,CR_cycle)
hold off
legend('eta charge','eta discharge','eta overall','CR')
title('Performance Indicators')