%% P=CPAR(Scenario)
% Create default P struct.
% Scenarios currently available: 
% - 'NYC2016'
% - 'NYC2018'
% 
% see also SCENARIOGENERATOR

function P=cpar(Scenario)

% scenario name (corresponding to associated file in 'data/scenarios/')
P.scenario=Scenario;    

% trip file associated with scenario
switch Scenario 
    case 'NYC2016'
        P.tripfile='NYC2016_Jan13-Mar16_10days';
    case 'NYC2018'
        P.tripfile='NY_trips_10wed_0103-0307_minutes';
end
P.scenarioid=1;     % which day of the trip file

P.eleproftype='TokyoDA-FY2017-Reduced'; % electricity price file in 'data/eleprices/'
P.eleprofseed=311;  % which day of the electricity file

% P.chargingstations=[4,4;5,6];

% simulation options
P.mpcpredict=1;
P.enlayeralg='aggregate';%'opti';
P.trlayeralg='simplified';%'opti';
P.v2g=true;
P.e=2;              % time step length in minutes
P.beta=15;

% technical parameters
P.m=10000;             % number of vehicles
P.battery=50;
P.chargekw=20;
P.consumption=0.15*30/60;   % consumption per minute: 0.15 kwh/km * 30 km/h / 60 min/h
P.cyclingcost=20000/2000;   % batterycost[yen/kWh] / lifelength[cycles] -> 10 yen/kWh

% operational parameters
P.initialsoc=0.7;
P.minsoc=0.2;
P.maxsoc=1;
P.v2gminsoc=0.5;
P.maxwait=10;

% % transport layer: opti
% Transport.thor=10;          % horizon (time steps)
% Transport.rho1=0.01;        % weight of secondary objective
% Transport.rho2=0.01;        % weight of charging objective for electricity price
% Transport.rho3=0.01;        % weight of charging objective for SOC
% Transport.rho4=0.000001;    % weight for fixed charge

% transport layer: simplified
Transport.tx=5;  % time steps
Transport.ts=15; % time steps
Transport.tr=15; % time steps

% energy layer: aggregate
Energy.mthor=48;         % macro horizon (macro time step)
Energy.bmin=0;
Energy.mminsoc=0.3;  % macro min soc
Energy.mmaxsoc=1;    % macro max soc

P.Energy=Energy;
P.Transport=Transport;


% % transport layer opti
% P.thor=10;          % horizon (time steps)
% P.rho1=0.01;        % weight of secondary objective
% P.rho2=0.01;        % weight of charging objective for electricity price
% P.rho3=0.01;        % weight of charging objective for SOC
% P.rho4=0.000001;    % weight for fixed charge
% 
% % energy layer: aggregate
% P.mthor=48;         % macro horizon (macro time step)
% P.tx=5;  % time steps
% P.ts=15; % time steps
% P.tr=15; % time steps
% P.bmin=0;
% P.mminsoc=0.3;  % macro min soc
% P.mmaxsoc=1;    % macro max soc

