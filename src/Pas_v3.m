
%%
% *This file is a draft and might very well contain mistakes. Please only use as inspiration.* 


%% INPUT Fix

% Mass flow:
M_tot = 750; % [kg/24h] --> [L/24h]
%M_tot = 2.4*35;
% variation mass flow to be here

% Piping:
pipe_D = 20; % [mm]
pipe_ins_thick = 100; % [mm]
pipe_ins_coef = 0.025; % [W/(K*m)] Kooltherm® Pipe Insulation 37 kg 

% Boudery conditions:
ambient_T = 273.15 + 24; % [K]
ambient_T = 273.15 + 30; % [K]
% ambient temp profile to be here

% Pathogen destruction/ contact chamber:
holding_T = 20; % [min]
holding_D = 20; % [mm]




%Massflow calc:
M_dot = M_tot/24/60/60; % [kg/s] --> [L/s] --> 0.001[m^3]
M_dot_m3 = M_dot*0.001; % [m^3/s]

%Holding calc:
h_t_s = holding_T *60; % [s] holding time in sec
holding_CA = (pi * (holding_D/2)^2)*10^-6 ; % [m^2], crossectional area
holding_L= h_t_s*M_dot_m3/holding_CA; % [m]


%% Input Variable

% System type: 0 = direct, 1 = water
type = 0;

%Gas Input:
V_gas = 1;  % [m^3/day] available gas per day 
% EC_gas = 0; % [J/m^3] Energy content gas per volume (pressure dependent?)
CH4_cont = 0.68; % methane content source: Hardeman thesis

%Burner Input:
eta_burner = 0.9; % [-] efficiency burner

% HX setting:
heat_coeff_1 = 300; % [W/K*m^2] source: VDI Wärmeatlas
heat_coeff_2 = 300; % [W/K*m^2] source: VDI Wärmeatlas
% or HX dep on flow see below

% Pasteurization setting:
pas_temp = 0;

%Type calc:
if type == 0
    HX1 = 'Pas_sim_v3/Type/Direct/HX1';
    Heat = 'Pas_sim_v3.Type.Direct.Heating.B.T';
elseif type == 1
    HX1 = 'Pas_sim_v3/Type/Water/HX1';
    HX2 = 'Pas_sim_v3/Type/Water/HX2';
    Heat = 'Pas_sim_v3.Type.Water.Heating.B.T';
end

%Energy input calc:
density_CH4 = 0.657; %[kg/m^3] @25C,1atm source: wikipedia
heat_comb_CH4 = 55.5; %[MJ/kg] 
EC_in_day = V_gas * density_CH4 * heat_comb_CH4 * CH4_cont * 1000000; % [J/day] energy in per day
EC_in = EC_in_day / 24 / 60 / 60; % [J/s] = [W] 
E_in = EC_in * eta_burner; % [W] energy after burner
KWH = density_CH4 * heat_comb_CH4 *1000/(60*60)*CH4_cont; % [kWh/m^3] (vs 6 source: Hardeman thesis)
%% HX calculations

% variable initialisation:

alpha = 0;              % Wärmeüberganszahl [W/(m²*K)] 
nu = 0;                 % Nusselt [-]
lambda_f = 0;           % Fluid heatconduction [W/(m*K)]
L_char = 0;             % characeristic lenght [m]

Pe = 0;                 % peclet [-]
Re = 0;                 % Reynolds [-]
Pr = 0;                 % Prandtl [-]

v = 0;                  % flow speed [m/s]
vis = 0;                % kin. viskosität [m^2/s]
a = 0;                  % temp leitzahl [m^2/s]

d = 0;                  % hydraulic diam 
d_ii = 0;               % inner diam inner pipe
d_oi = 0;               % outer diam inner pipe
d_io = 0;               % inner diam outer pipe
L_s = 0;                % length hx section
K = 0;                  % K factor
Pr_f = 0;               % Prantl fluid
Pr_w = 0;               % Prantl wand



% calc nu for laminar flow:

nu = (49.37 + (1.615*(Pe*d/h)^(1/3)-0.7)^3)^(1/3)*K;

%% HX v2 Forbis

%TL1 shell, TL2 tube

%Input variables tl1:
TL1_L = 0.762;                % length of flow path for tl 1 [m]    source: Forbis hx
D1_i = 125;                   % inner diameter tl1 [mm]             source: Forbis hx
D1_o = 127;                   % outer diameter tl1 [mm]             source: Forbis hx
fouling_factor_1 = 1.41       % [K*m^2/kW]                          source: VDI Wärmeatlas

%Input variables tl2:
TL2_L = TL1_L;                  % length of flow path for tl2 [m]   source: Forbis hx
D2_i = 55;                      % inner diameter tl2 [mm]           source: Forbis hx
D2_o = 57;                      % outer diameter tl2 [mm]           source: Forbis hx
fouling_factor_2 = 1.41         % [K*m^2/kW]                        source: VDI Wärmeatlas

%calc_temp variables:
D1_im = D1_i * 0.001;
D1_om = D1_o * 0.001;
D2_im = D2_i * 0.001;
D2_om = D2_o * 0.001;
R1_i = D1_im/2;
R1_o = D1_om/2;
R2_i = D2_im/2;
R2_o = D2_om/2;

%Output variables TL1:
min_flow_area_1 = (R1_i^2 - R2_o^2)*pi;             % [m^2]
Dh_press_1 = D1_im - D2_om;                         % [m]
liquid_volume_1 = min_flow_area_1 * TL1_L;          % [m^3]
heat_transfer_area_1 = 2*pi*R2_o*TL1_L;             % [m^2]

%Output variables TL2:
min_flow_area_2 = R2_i^2 * pi; 
Dh_press_2 = D2_im;
liquid_volume_2 = min_flow_area_2 * TL2_L;
heat_transfer_area_2 = 2*pi*R2_i*TL2_L;
%% HX v2 ETH Zürich

%TL1 shell, TL2 tube

%Input variables tl1:
TL1_L = 4;                      % length of flow path for tl 1 [m]    source: ETH Zürich hx
D1_i = 48;                      % inner diameter tl1 [mm]             source: ETH Zürich hx
D1_o = 50;                      % outer diameter tl1 [mm]             source: ETH Zürich hx
fouling_factor_1 = 1.41;        % [K*m^2/kW]                          source: VDI Wärmeatlas

%Input variables tl2:
TL2_L = TL1_L;                  % length of flow path for tl2 [m]   source: ETH Zürich hx
D2_i = 18;                      % inner diameter tl2 [mm]           source: ETH Zürich hx
D2_o = 20;                      % outer diameter tl2 [mm]           source: ETH Zürich hx
fouling_factor_2 = 1.41;        % [K*m^2/kW]                        source: VDI Wärmeatlas

%calc_temp variables:
D1_im = D1_i * 0.001;
D1_om = D1_o * 0.001;
D2_im = D2_i * 0.001;
D2_om = D2_o * 0.001;
R1_i = D1_im/2;
R1_o = D1_om/2;
R2_i = D2_im/2;
R2_o = D2_om/2;

%Output variables TL1:
min_flow_area_1 = (R1_i^2 - R2_o^2)*pi;             % [m^2]
Dh_press_1 = D1_im - D2_om;                         % [m]
liquid_volume_1 = min_flow_area_1 * TL1_L;          % [m^3]
heat_transfer_area_1 = 2*pi*R2_o*TL1_L;             % [m^2]

%Output variables TL2:
min_flow_area_2 = R2_i^2 * pi; 
Dh_press_2 = D2_im;
liquid_volume_2 = min_flow_area_2 * TL2_L;
heat_transfer_area_2 = 2*pi*R2_i*TL2_L;

%% HX v2 Plate Gasketed Plate Heat Exchanger T2-BFG 35 plates ISO 7 - R 3/4 connection

%TL cold , TL2 warm

%Input variables tl1:
TL1_L = 0.39*17+0.23;                      % length of flow path for tl 1 [m]    source: approximation
fouling_factor_1 = 1.41;        % [K*m^2/kW]                          source: VDI Wärmeatlas

%Input variables tl2:
TL2_L = TL1_L;                  % length of flow path for tl2 [m]   source: ETH Zürich hx
fouling_factor_2 = 1.41;        % [K*m^2/kW]                        source: VDI Wärmeatlas

Plate_con_d = 0.02; %diameter conection [m]

%Output variables TL1:
min_flow_area_1 = Plate_con_d^2*pi;             % [m^2]
Dh_press_1 = 0.00472;                         % [m] source: chat gpt and data sheet
liquid_volume_1 = min_flow_area_1 * TL1_L;          % [m^3]
heat_transfer_area_1 = 1.26;             % [m^2]

%Output variables TL2:
min_flow_area_2 = Plate_con_d^2*pi; 
Dh_press_2 = 0.00472;
liquid_volume_2 = min_flow_area_2 * TL2_L;
heat_transfer_area_2 = 1.26;


%% Simulation command

open_system('Pas_sim_v3.slx');
E_in_range = 100 : 50 : 300;
range_length = length(E_in_range);

for i = range_length:-1:1
    in(i) = Simulink.SimulationInput('Pas_sim_v3');

 % Gas / Energy in
 % in = in.setBlockParameter('Pas_sim_v3/Q_in',...
      % 'value',num2str(E_in));
    in(i) = in(i).setBlockParameter('Pas_sim_v3/Q_in',...
       'value',num2str(E_in_range(i)));

 % %HX
 %     in(i) = in(i).setBlockParameter(HX1, ...
 %     'heat_coeff_1', num2str(heat_coeff_1));
 %     in(i) = in(i).setBlockParameter(HX1, ...
 %     'heat_coeff_2', num2str(heat_coeff_2));
 %     if type == 1
 %        in(i) = in(i).setBlockParameter(HX2, ...
 %         'heat_coeff_1', num2str(heat_coeff_1));
 %         in(i) = in(i).setBlockParameter(HX2, ...
 %         'heat_coeff_2', num2str(heat_coeff_2));
 %     end
end
 out = parsim(in,'ShowSimulationManager','on','ShowProgress','on');

%% Data analysis
pas_temp_matrix = zeros(range_length,2)

for a = 1 : 1 : range_length % FINISH THIS LOOP
    H_exit_t = getElement(out(1,a).xout,'Pas_sim_v3.Type.Water_Direct1.HX2.thermal_liquid_2.T_A').Values;
    H_exit_t_table = timeseries2timetable(H_exit_t);
    H_exit_t_table_sec = retime(H_exit_t_table,"secondly", "linear");
    H_exit_t_lenght = height(H_exit_t_table_sec);

    for i = 2001:2000:H_exit_t_lenght
        if ismembertol(H_exit_t_table_sec{i,:},H_exit_t_table_sec{i-2000,:},0.0000001)
            pas_temp = H_exit_t_table_sec{i,:}
        end
    end
    pas_temp_matrix(a,1) = E_in_range(a);
    pas_temp_matrix(a,2) = pas_temp;
end

%%
figure
plot(pas_temp_matrix(:,1), pas_temp_matrix(:,2))