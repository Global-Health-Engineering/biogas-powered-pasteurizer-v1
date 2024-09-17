%% input parameters
d_s = 0.127; %[m]
d_to = 0.057+0.0032; %[m]
d_ti = 0.057; %[m]
L = 0.762; %[m]
W_p = 2.4*10^(-3); %[L/p] wastewater per person
P= 26; %[people]
n_u = 0.9131*10^(-6); %[m^2/s] @24C
a = 0.1450*10^(-6); %[m^2/s] @24C
lam = 604.87*10^(-3); %[W/mK]
lam_w = 50 ; %[W/mK]
Nu_t = 3.66; % [-]

%% calations:

V_dot = (W_p*P) / 86400;

Pr = n_u / a;

w_s = V_dot / (pi* ((d_s/2)^2-(d_to/2)^2) );

d_h = d_s - d_to;

Re = w_s * d_h / n_u;

Nu_1s = 3.66 +1.2*(d_s/d_to)^(-0.8);

Nu_2s = 1.615*(1+0.14*(d_to/d_s)^(-0.5))*(Re*Pr*d_h/L)^(1/3);

Nu_3s = ( (2/(1+22*Pr))^(1/6) ) * (Re*Pr*d_h/L)^0.5;

Nu_s = (Nu_1s^3 + Nu_2s^3 +Nu_3s^3)^(1/3);

alpha_s = Nu_s * lam / d_h;

alpha_t = Nu_t * lam / d_ti;

k = 1 / ( (d_to/(alpha_s*d_s)) + (d_ti * log(d_ti/d_s)/(2*lam_w) + (1/alpha_t) )) ;

area = pi*d_to*L;

a_v = pi*d_to*L / (V_dot);
