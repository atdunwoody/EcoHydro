%% TOY MODEL

function [Ea,QF,R,QS,QT,Sf,Su,Ss,St,AL,IE,SE]= toymodel(DATA,PAR)

[M,~]   = size(DATA);
P       = DATA(:,1); 
Ep      = DATA(:,2);
PAR(5) = 1;   % split between R and Qfs [0-1];


    mir        = PAR(1);  % Maximum infiltration rate 
    Su_max     = PAR(2);  % Unsaturated zone storage capacity [50-300]
    Ts         = PAR(3);  % time parameter for slowflow [20-100]
    Tf         = PAR(4);  % time parameter for quickflow [1-3]
    beta       = PAR(5);  % Split between Recharge and overland flow [SET IT TO = 1]

% Initialization
QF=zeros(M,1);
QS=zeros(M,1);
SE=zeros(M,1);

QT=zeros(M,1);
R=zeros(M,1);
Sf=zeros(M,1);
Su=zeros(M,1);
Ss=zeros(M,1);
Ea=zeros(M,1);
KGE=zeros(M,1);
AL=zeros(M,1);
IE = zeros(M,1);
SE = zeros(M,1);

% Initial values
S0    = 0.2; % initial condition
Su_dt = S0;
Ss_dt = S0;
Sf_dt = S0;

% mir = 0.2;

%% MAIN FLOW ROUTINE
    
for t=1:M
             if P(t)>0
             AL(t) = 1 - (1 - exp(-P(t)/mir) )/(P(t)/mir);
             else
             AL(t) = 0;
             end
             
            %Unsaturated zone water balance    
            [r,se,E,Su_dt] = SU_eq(Su_dt,P(t),Ep(t),Su_max,AL(t),beta);
            
            R(t)  = r;       % recharge
            Ea(t) = E;      % evaporation 
            Su(t) = Su_dt;
            SE(t) = se;     % Saturation Excess Overland Flow Input [NOT CONSIDERED WHEN beta = 1]
            
            
            
            %Saturated zone
            [Qs,Ss_dt] = SS_eq(Ss_dt,r,Ts);
            
            QS(t) = Qs;     % slowflow
            Ss(t) = Ss_dt;  

end


for t=1:M
    
    %Unsaturated zone    
    [Qf,Sf_dt,ie]   = SF_eq(Sf_dt,P(t),AL(t),Tf,SE(t)); 

    QF(t) = Qf;         % quickflow
    Sf(t) = Sf_dt;      
    IE(t) = ie;

end

% model's outputs
Ea = Ea;                
QS = QS;
QF = QF;
R = R;
QT = QF + QS;
Sf = Sf;
Su = Su;
Ss = Ss;
AL = AL;
 
St                               = Su + Ss;

end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fast Flow Reservoir (Sf)

function [Qf,S_dt,ie]= SF_eq(S0,P,alpha,Tf,SE)
 
 ie    = (P*alpha);
 S_dt  = S0 + (P*alpha) + SE - (S0/Tf);
 Qf    = S0/Tf;
end

%% Unsaturated Zone Reservoir (Su)


function [r,se,E,S_dt]= SU_eq(Su_0,P,Ep,Sumax,alpha,beta)

 % NEW - 04-12-2022
 Su_0 = Su_0 + (P*(1-alpha));
 
 
 if Su_0<Sumax
      R = 0;
 else
      R = Su_0 - Sumax;
 end

 Su_0 = Su_0 - R;
 
 E = (Ep*(Su_0/Sumax));
 if E>Ep
    E=Ep;
 end
 Su_0 = Su_0 - E;
 
 S_dt = Su_0; 
 r  = R*beta;
 se = R*(1-beta);

end
%% Saturated Zone Reservoir (Ss)

function [Qs,S_dt]= SS_eq(S0,R,Ts)
 
 S_dt = S0 + R - (S0/Ts);
 Qs   = S0/Ts;
end




