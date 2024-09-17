%% Read Me:
% This file can be used togheter with pas_sim_v6 to calculate the results
% from the global simulations for the holding analysis in the thesis.

%% Input manually Variable
%This and all section bellow need to be run 4 times, for the four different
%holding times. 


%changed from 30, 60, 120 and 180 min
holding_T = 120; % [min] Holding retention time, source 120 min: Presentation: Three years of field experience piloting the Anaerobic Digestion Pasteurization Latrine by Forbis-Stokes et al.




%% INPUT Fix
%These parameters were fixed and not varied in the global simulations

% Mass/Volume flow:
M_tot = 750; %  [L/24h] based on: thesis Hardeman https://doi.org/10.3929/ethz-b-000587768



% Piping:
pipe_D = 20; % [mm]
pipe_ins_thick = 100; % [mm] twice the insulation thickness because of diameter
pipe_ins_coef = 0.025; % [W/(K*m)] Kooltherm® Pipe Insulation 37 kg https://www.kingspan.com/content/dam/kingspan/kti/products/general-gb-and-ireland/kingspan-kooltherm-pipe-insulation-data-sheet-en-ie-gb.pdf
pipe_L = 0.5; % [m]

% Boundary conditions:
ambient_T = 273.15 + 24; % [K] ambient temperature


% Pathogen inactivation/ contact chamber:
holding_D = 40; % [mm] Diamiter of holding pipe

% HX coefficients:
heat_coeff_1 = 100; % [W/K*m^2] source:  own calculations or VDI Wärmeatlas https://doi.org/10.1007/978-3-662-52989-8
heat_coeff_2 = 100; % [W/K*m^2] source: own calculations or VDI Wärmeatlas https://doi.org/10.1007/978-3-662-52989-8


% Hx calc:
heat_transfer_area_1 = 0.144112394842532* (750/(2.4*35)); %surface area from ADPL corrected for different massflow https://doi.org/10.1089/ees.2016.0148
heat_transfer_area_2 = 0.144112394842532* (750/(2.4*35)); %symmetrical


%Massflow calc:
M_dot = M_tot/24/60/60; % [L/s] --> 0.001[m^3]
M_dot_m3 = M_dot*0.001; % [m^3/s]

%Holding calc:
h_t_s = holding_T *60; % [s] holding time in sec
holding_CA = (pi * (holding_D/2)^2)*10^-6 ; % [m^2], cross-sectional area
holding_L= h_t_s*M_dot_m3/holding_CA; % [m] length of holding pipe



%% Automatic Variable 
%This section defines the range of energy inputs over wich was varied in
%the global simulations

%Gas Input:v
V_gas1 = 1.5;  % [m^3/day] available gas per day, lower limit
V_gas2 = 5;    % [m^3/day] available gas per day, upper limit
CH4_cont = 0.68; % methane content source: thesis Hardeman https://doi.org/10.3929/ethz-b-000587768


%Energy input calc:
density_CH4 = 0.657; %[kg/m^3] @25C,1atm source: wikipedia https://en.wikipedia.org/wiki/Methane
heat_comb_CH4 = 55.5; %[MJ/kg] Heat of combustion source: ??
EC_in_day1 = V_gas1 * density_CH4 * heat_comb_CH4 * CH4_cont * 1000000; % [J/day] energy in per day
EC_in_day2 = V_gas2 * density_CH4 * heat_comb_CH4 * CH4_cont * 1000000; % [J/day] energy in per day
EC_in1 = EC_in_day1 / 24 / 60 / 60; % [J/s] = [W] 
EC_in2 = EC_in_day2 / 24 / 60 / 60; % [J/s] = [W] 
E_in1 = EC_in1 * 0.55; % [W] energy after burner, lower limit. Heater efficiency of 0.55
E_in2 = EC_in2 * 0.86; % [W] energy after burner, upper limit. Heater efficiency of 0.86
KWH = density_CH4 * heat_comb_CH4 *1000/(60*60)*CH4_cont; % [kWh/m^3] (vs 6 source: hardman thesis)

%% HX variables for pressure calculations, based on ADPL heat exchanger, dont have an influence on heat exchange
%These variables need to be defined for simulink to run, but should not
%effect pasteurisation temperature.


%TL1 shell, TL2 tube

%Input variables tl1:
TL1_L = 0.762;                % length of flow path for tl 1 [m]    source: ADPL heat exchanger https://doi.org/10.1089/ees.2016.0148
D1_i = 127;                   % inner diameter tl1 [mm]             source: ADPL heat exchanger https://doi.org/10.1089/ees.2016.0148
D1_o = 127 + 3.2;             % outer diameter tl1 [mm]             source: ADPL heat exchanger https://doi.org/10.1089/ees.2016.0148
fouling_factor_1 = 1.41   ;    % [K*m^2/kW]                         source: VDI Wärmeatlas https://doi.org/10.1007/978-3-662-52989-8

%Input variables tl2:
TL2_L = TL1_L;                  % length of flow path for tl2 [m]   source: ADPL heat exchanger https://doi.org/10.1089/ees.2016.0148
D2_i = 57;                      % inner diameter tl2 [mm]           source: ADPL heat exchanger https://doi.org/10.1089/ees.2016.0148
D2_o = 57 + 3.2;                % outer diameter tl2 [mm]           source: ADPL heat exchanger https://doi.org/10.1089/ees.2016.0148
fouling_factor_2 = 1.41  ;       % [K*m^2/kW]                       source: VDI Wärmeatlas https://doi.org/10.1007/978-3-662-52989-8

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
%heat_transfer_area_1 = 2*pi*R2_o*TL1_L;             % [m^2] this variable was used to define the variable heat_transfer_area_1 and heat_transfer_area_2 in "input fix" 

%Output variables TL2:
min_flow_area_2 = R2_i^2 * pi; 
Dh_press_2 = D2_im;
liquid_volume_2 = min_flow_area_2 * TL2_L;
%heat_transfer_area_2 = 2*pi*R2_i*TL2_L; this variable was not used to define the variable heat_transfer_area_1 and heat_transfer_area_2 in "input fix" to keep the heat exhanger symetric. However, the change would be negligible.


%% Simulation command
% Section needs to be run three times

open_system('Pas_sim_v6.slx');
E_in_range = E_in1 : 50 : E_in2; %creates vertor from the lower to the upper energy input
range_length = length(E_in_range);

for i = range_length:-1:1
    in(i) = Simulink.SimulationInput('Pas_sim_v6');

%Defines multiple simulations with the different energy inputs:
    in(i) = in(i).setBlockParameter('Pas_sim_v6/Q_in',...
       'value',num2str(E_in_range(i)));
    %in(i) = in(i).setBlockParameter('Pas_sim_v6/Type/Direct/Holding',...
       %'pipe_lenght',num2str(E_in_range(i)));
end
 out = parsim(in,'ShowSimulationManager','on','ShowProgress','on',TransferBaseWorkspaceVariables='on'); %Runs the before defined simulations and returns results.

%% Data analysis
% Section needs to be run three times


pas_temp_matrix4 = zeros(range_length,2); %first run, name this variable 'pas_temp_matrix1', second run: 'pas_temp_matrix2', third run: 'pas_temp_matrix3'
t90_matrix4 = zeros(range_length,2); %first run, name this variable '...1', second run: '...2', third run: '...3'
t90=0;

for a = 1 : 1 : range_length 
    H_exit_t = getElement(out(1,a).xout,'Pas_sim_v6.Type.Direct.Holding.B.T').Values; %reads temp value after holding
    H_exit_t_table = timeseries2timetable(H_exit_t); %changes matlab storage type
    H_exit_t_table_sec = retime(H_exit_t_table,"secondly", "linear"); % linearly interpolates the resulst from the simulation with a clearly defined time interval
    H_exit_t_lenght = height(H_exit_t_table_sec);
    
    for i = 2001:2000:H_exit_t_lenght
        if ismembertol(H_exit_t_table_sec{i,:},H_exit_t_table_sec{i-2000,:},0.00001) %checks if simulation is stationary 
            pas_temp = H_exit_t_table_sec{i,:};
        end
    end
    for ii = 21:20:H_exit_t_lenght
        if ismembertol(H_exit_t_table_sec{ii,:}-273.15,(pas_temp-273.15)*0.9,0.001) %checks if simulation arrived at 0.9 of stationary temp
            t90 = ii;
        end
    end

    pas_temp_matrix4(a,1) = E_in_range(a); %first run, change variable '...1', second run: '...2', third run: '...3'
    pas_temp_matrix4(a,2) = pas_temp - 273.15; % store final pasteurisation temperatures in degree celcius. %first run, change variable '...1', second run: '...2', third run: '...3'
    t90_matrix4(a,1) = E_in_range(a); %first run, change variable '...1', second run: '...2', third run: '...3'
    t90_matrix4(a,2) = t90; %first run, change variable '...1', second run: '...2', third run: '...3'

end

%% Figure Pasteurisation temp vs Energy
% Section needs to be run 1 time


figure
plot(pas_temp_matrix1(:,1), pas_temp_matrix1(:,2)) %takes the results of third run
title('1-1')
xlabel('Energy in [W]', 'Interpreter', 'tex');
ylabel('Temperature [°C]', 'Interpreter', 'tex');
hold on
plot(pas_temp_matrix2(:,1), pas_temp_matrix2(:,2)) %takes the results of second run
plot(pas_temp_matrix3(:,1), pas_temp_matrix3(:,2)) %takes the results of first run
plot(pas_temp_matrix4(:,1), pas_temp_matrix4(:,2)) %takes the results of first run
hold off


%% Figure T90 time vs Energy
% Section needs to be run 1 time


figure
plot(t90_matrix1(:,1), t90_matrix1(:,2)) %takes the results of third run
title('t90 vs energy in')
xlabel('Energy in [W]', 'Interpreter', 'tex');
ylabel('Time [s]', 'Interpreter', 'tex');
hold on
plot(t90_matrix2(:,1), t90_matrix2(:,2)) %takes the results of second run
plot(t90_matrix3(:,1), t90_matrix3(:,2)) %takes the results of first run
plot(t90_matrix4(:,1), t90_matrix4(:,2)) %takes the results of first run
hold off