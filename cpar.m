%% P=CPAR(Scenario)
% Create default P struct.
% Scenarios currently available: 
% - 'NYC2016'
% - 'NYC2018'
% - 'NYC2016-small'
% - 'NYC2016-small2'
% 
% see also SCENARIOGENERATOR

function P=cpar(Scenario,trlayer)

% scenario name (corresponding to associated file in 'data/scenarios/')
P.scenario=Scenario;

if nargin<2
    trlayer='simplified';
end

% trip file associated with scenario
switch Scenario 
    case 'NYC2016'
        P.tripfile='NYC2016_Jan13-Mar16_10days';
    case 'NYC2018'
        P.tripfile='NY_trips_10wed_0103-0307_minutes';
    case 'NYC2016-small'
        P.tripfile='NYC2016-small_13Jan';
    case 'NYC2016-small2'
        P.scenario='NYC2016-small';
        P.tripfile='NYC2016-small2_13Jan';
end
P.scenarioid=1;     % which day of the trip file

P.eleproftype='TokyoDA-FY2017-Reduced'; % electricity price file in 'data/eleprices/'
P.eleprofseed=311;  % which day of the electricity file

% P.chargingstations=[4,4;5,6];

% simulation options
P.mpcpredict=true;
P.modechoice=false;
P.v2g=true;
P.e=2;              % time step length in minutes
P.beta=15;
P.m=10000;             % number of vehicles

% technical parameters
Tech.battery=50;
Tech.chargekw=20;
Tech.consumption=0.15*30/60;   % consumption per minute: 0.15 kwh/km * 30 km/h / 60 min/h
Tech.cyclingcost=20000/2000;   % batterycost[yen/kWh] / lifelength[cycles] -> 10 yen/kWh

% operational parameters
Operations.initialsoc=0.7;
Operations.minsoc=0.2;
Operations.maxsoc=1;
Operations.v2gminsoc=0.5;
Operations.maxwait=10;

P.trlayeralg=trlayer;

% % transport layer: opti
if strcmp(P.trlayeralg,'opti')
    Transport.thor=8;          % horizon (time steps)
    Transport.rho1=0.01;        % weight of secondary objective
    Transport.rho2=0.01;        % weight of charging objective for electricity price
    Transport.rho3=0.01;        % weight of charging objective for SOC
    Transport.rho4=0.000001;    % weight for fixed charge
end

% transport layer: simplified
if strcmp(P.trlayeralg,'simplified')
    Transport.tx=5;  % time steps
    Transport.ts=15; % time steps
    Transport.tr=15; % time steps
    Transport.bmin=0;
end

% energy layer: aggregate
P.enlayeralg='aggregate';
Energy.mthor=48;         % macro horizon (macro time step)
Energy.mminsoc=0.3;  % macro min soc
Energy.mmaxsoc=1;    % macro max soc

% append to main struct
P.Tech=Tech;
P.Operations=Operations;
P.EnergyLayer=Energy;
P.TransportLayer=Transport;

% case-dependent adjustments
switch Scenario 
    case 'NYC2016-small'
        P.m=30; 
    case 'NYC2016-small2'
        P.m=5;
end

