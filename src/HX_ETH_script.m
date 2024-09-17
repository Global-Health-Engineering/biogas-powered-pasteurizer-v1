%% input parameters
d_s = 0.048; %[m]
d_to = 0.02; %[m
d_ti = 0.018; %[m]
L = 4; %[m]
V = 0.750; %[m^3/d]
n_u = 0.9131*10^(-6); %[m^2/s] @24C kinimatic viscocity
a = 0.1450*10^(-6); %[m^2/s] @24C
lam = 604.87*10^(-3); %[W/mK]
lam_w = 50 ; %[W/mK]
Nu_t = 3.66; % [-]

%% calations:

V_dot = V / 86400; %volume flow

Pr = n_u / a; %prantel number

w_s = V_dot / (pi* ((d_s/2)^2-(d_to/2)^2) ); %flow speed

d_h = d_s - d_to; % hydraulic diamiter

Re = w_s * d_h / n_u; % renolds number

Nu_1s = 3.66 +1.2*(d_s/d_to)^(-0.8);

Nu_2s = 1.615*(1+0.14*(d_to/d_s)^(-0.5))*(Re*Pr*d_h/L)^(1/3);

Nu_3s = ( (2/(1+22*Pr))^(1/6) ) * (Re*Pr*d_h/L)^0.5;

Nu_s = (Nu_1s^3 + Nu_2s^3 +Nu_3s^3)^(1/3);

alpha_s = Nu_s * lam / d_h;

alpha_t = Nu_t * lam / d_ti;

k = 1 / ( (d_to/(alpha_s*d_s)) + (d_ti * log(d_ti/d_s)/(2*lam_w) + (1/alpha_t))) ;

area = pi*d_to*L;

R_w = (d_to - d_ti)/ (lam_w*area);

a_v = pi*d_to*L / (V_dot);

%% flow to achieve turbulance 

w_s_critical = 2300 * n_u / d_h;
V_dot_critical = w_s_critical * (pi* ((d_s/2)^2-(d_to/2)^2) );
V_critical = V_dot_critical * 86400; % 9.6908 m^3
