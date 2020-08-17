%% P=CPAR(Scenario[,trlayer])
% Create default P struct.
% Scenarios currently available: 
% - 'NYC2016'
% - 'NYC2018'
% - 'NYC2016-small'
% - 'NYC2016-small2'
% - 'Tokyo189'
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
        % P.tripfile='NYC2016_Jan13-Mar16_10days';
        % P.tripfolder='NYC2016';
        P.tripfolder='NYC2016-nodes';
    case 'NYC2018'
        % P.tripfolder='NYC2018_10wed';
        % P.tripfile='NY_trips_10wed_0103-0307_minutes';
        P.tripfolder='NYC2018_10wed-nodes';
    case 'NYC2016-small'
        P.tripfile='NYC2016-small_13Jan';
    case 'NYC2016-small2'
        P.scenario='NYC2016-small';
        P.tripfile='NYC2016-small2_13Jan';
    case 'Tokyo189'
        P.tripfile='Tokyo2008_1day_48k';
    case 'Munich'
        P.tripfile='Munich';
end
P.tripday=1;     % which day of the trip file

P.gridfile='TokyoDA-FY2017-Reduced'; % electricity price file in 'data/eleprices/' ($/MWh; half-hourly, 48x365)
P.gridday=311;  % which day of the electricity file

% P.chargingstations=[4,4;5,6];

% simulation options
P.mpcpredict=true;      % perfect prediction?
P.modechoice=false;     % mode choice?
P.e=2;                  % time step length in minutes
P.beta=30;              % frequency of energy layer simulation (in minutes)
P.m=10000;              % number of vehicles
P.carbonprice=0;        % $/ton

% technical parameters
Tech.battery=50;                % kWh
Tech.chargekw=20;               % kW
Tech.consumption=0.15*30/60;    % kWh/min   (0.15 kwh/km * 30 km/h / 60 min/h)
Tech.cyclingcost=100/4000;%200/2000;      % $/kWh     (batterycost[$/kWh] / lifelength[cycles])

% operational parameters
Operations.initialsoc=0.7;
Operations.minsoc=0.2;
Operations.maxsoc=1;
Operations.v2g=true;
Operations.v2gminsoc=0.5;
Operations.maxwait=10; % minutes

P.trlayeralg=trlayer;

% % % transport layer: opti
% if strcmp(P.trlayeralg,'opti')
%     Transport.thor=8;           % horizon (time steps)
%     Transport.rho1=0.01;        % weight of secondary objective
%     Transport.rho2=0.01;        % weight of charging objective for electricity price
%     Transport.rho3=0.01;        % weight of charging objective for SOC
%     Transport.rho4=0.000001;    % weight for fixed charge
% end

% transport layer: simplified
if strcmp(P.trlayeralg,'simplified')
    Transport.tx=10;  % minutes
    Transport.ts=24; 
    Transport.tr=20; 
    Transport.bmin=0;
end

% energy layer: aggregate
P.enlayeralg='aggregate';
Energy.mthor=1440;      % macro horizon (minutes)
Energy.extrasoc=0.25;  % extra soc for energy layer to account for aggregate uncertainty

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
    case 'Tokyo189'
        P.m=3500;
end

