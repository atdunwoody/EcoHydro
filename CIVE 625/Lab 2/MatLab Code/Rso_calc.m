
function [Ra, Rso] = Rso_calc(Dates,Lat,Lon,Elev)



%% Make Julian dates
years = unique(Dates(:,1));
jday=0;
for i=1:length(years)
    
    aux      = find(Dates(:,1) == years(i));
    jday_new = [1:1:length(aux)]';
    jday     = [jday;jday_new];

end
jday(1) = [];

%% Calculate Ra or Sod Incoming radiation without atmosphere
clear Ra

dr      = 1+0.033.*cos(2.*pi.*jday/365);
delta   = 0.409.*sin(2.*pi.*jday/365-1.39);
Lat     = Lat.*pi/180;
Long    = Lon.*pi/180;
Z       = Elev;



for i=1:length(Lat)
    
    ws                = (acos(-tan(delta).*tan(Lat(i))));

%   Ra(:,i)           = 24*60/pi*0.0820.*dr.*(ws.*(sin(Lat(i)).*sin(delta))+sin(ws).*(cos(Lat(i)).*cos(delta))); %MJ. m-2 .day-1
    Ra(:,i)           = 24*60/pi*(0.08202).*dr.*(ws.*(sin(Lat(i)).*sin(delta))+sin(ws).*(cos(Lat(i)).*cos(delta))); %MJ. m-2 .day-1

    Rso(:,i)          = (0.75+2*10^-5*Z(i)).*Ra(:,i) ;                                                            %% Calculate Rso, clear sky radiation
end
