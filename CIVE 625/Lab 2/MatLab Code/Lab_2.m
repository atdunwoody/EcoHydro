clear all
clc 
close all

cd 'C:\Users\ameira\OneDrive - Colostate\4 Teaching\CIVE 625 - Ecohydrology\Labs\Lab 2'
%% 1 - Load Data
Q              = csvread('Q_mm.csv');   % Q in mm 
Dates          = Q(:,1:3); 
Precip         = csvread('Precip.csv'); % P in mm 
e              = csvread('e.csv');      % e in kPa
u2             = csvread('u2.csv');     % Wind Speed at 2m in m/s
S_in           = csvread('S_in.csv');   % Solar radiation in W/m2
Temp           = csvread('Temp.csv');   % Temp in Celsius
Lat_Lon_A_Z    = csvread('Lat_Lon_Area_Z.csv');   % Temp in Celsius

Q(:,1:3)       = [];
Precip(:,1:3)  = [];
e(:,1:3)       = [];
u2(:,1:3)      = [];
S_in(:,1:3)    = [];
Temp(:,1:3)    = [];

%% 2.1 - Generate Daily PET: Calculate Net Longwave

[Ra,Rso]    = Rso_calc(Dates,Lat_Lon_A_Z(:,1),Lat_Lon_A_Z(:,2),Lat_Lon_A_Z(:,4));    % Rso = Clear sky solar radiation MJ/(m2.day);

Rso        = Rso./(10^-6*60*60*24); % Convert to W/m2
Ra         = Ra./(10^-6*60*60*24); % Convert to W/m2

from  = 100; % start plot at
to    = 300; % end plot at
chose = 15;  % pick catchment

figure(1)
clf(1)
plot(Ra(from:to,chose),'.b'); hold on
plot(Rso(from:to,chose),'.k'); hold on
plot(S_in(from:to,chose),'.-r'); hold on
legend('Extraterrestrial solar radiation','Clear sky solar radiation','S in')
ylabel('W/m^2')
xlabel('days')

f        = S_in./Rso;
e_prime  = 0.34 - 0.14*(e./1000).^0.5;
L_net    = -f.*e_prime.*5.67*10^-8.*(Temp+273.15).^4; % Equation 5.22 Terrestrial Hydrometeorology (TH)


%% 2.2 - Generate Daily PET:
cp       = 1.1013;                                        % specific heat at constant pressure for air kJ/(kg.K)
lambda   = 10^3*(2.501 - 0.002361.*Temp);                 % Latent Heat of Vaporization (kJ/kg)
rho_air  = 1.23;
ra       = 208./u2;
rs       = 70; 
Rn_RC    = S_in.*(1-0.23) + L_net;                        % Net Radiation over a Reference Crop, in W/m2
e_sat    = 0.6108.*exp( (17.27.*Temp)./(237.3 + Temp) );  % Saturated Vapor Pressure in kPa
D        = e_sat - e;                                     % Vapor Pressure Deficit in kPa
Delta    = 4098.*e_sat./(237.3 + Temp).^2;                % Slope of esat versus temp curve (kPa/C)
Z        = Lat_Lon_A_Z(:,4);                              % Elevation (m)
Press    = 101.3*((293-0.0065.*Z)/293).^5.26;             % Pressure as function of elevation (kPa)
gamma    = cp.*Press'./(0.622.*lambda);                   % Psychrometric Constant (kPa/degC)

Wm2_mm  = (lambda.*10^3).^-1.*86400;                      % Conversion Factor from W/m2 to mm/day

E_RC    = ( Delta.*Rn_RC + rho_air.*cp.*D./ra )./( Delta + gamma.*( 1 + rs./ra) );

E_RC_mm  = E_RC.*Wm2_mm;
Rn_RC_mm = Rn_RC.*Wm2_mm;

%% Plot Example
from  = 500; % start plot at
to    = 800; % end plot at
chose = 13;  % pick catchment

figure(1)
clf(1)
subplot(1,2,1)
plot(E_RC(from:to,chose),'.b'); hold on
plot(Rn_RC(from:to,chose),'.-k'); hold on
legend('Ref. Crop ET','Net Radiation')
ylabel('W/m2')
xlabel('days')

subplot(1,2,2)
plot(E_RC_mm(from:to,chose),'.b'); hold on
plot(Rn_RC_mm(from:to,chose),'.-k'); hold on
legend('Ref. Crop ET','Net Radiation')
ylabel('mm')
xlabel('days')


%% Run Hydrologic Model
pick_catchment = 8;

INPUT(:,1)  = Precip(:,pick_catchment);
INPUT(:,2)  = E_RC_mm(:,pick_catchment);

PAR(1) = 50;  % % Maximum infiltration rate  [300]
PAR(2) = 450;   % Su_max - total water capacity [50-300]
PAR(3) = 5;    % Ts time - parameter for slowflow [20-100]
PAR(4) = 1;     % Tf time = parameter for quickflow [1-3]

[Ea,QF,R,QS,QT,Sf,Su,Ss,St,AL,IE,SE] = toymodel(INPUT,PAR); 
%%
Su_max   = PAR(2);

figure(1)
clf(1)
subplot(1,2,1)
plot(Q(from:to,pick_catchment),'-b'); hold on
plot(QT(from:to),'-k'); hold on
legend('Observed','Simulated')
ylabel('mm')
xlabel('days')
% 
subplot(1,2,2)
plot(Su(from:to)./Su_max,'.-b'); hold on
ylabel('[-]')

yyaxis right
plot(Ea(from:to),'.-r'); hold on

legend('Su/Su_{max}','ET')
ylabel('mm')
xlabel('days')
% 





