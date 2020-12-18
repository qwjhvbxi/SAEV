% scenario
P.scenario='NYC2016';
P.tripfolder='NYC2016-nodes';
P.tripday=13;     % which day of the trip file

P.gridfile='NY_DA_2016'; % electricity price file in 'eleprices/' ($/MWh; half-hourly, 48x365)
P.gridday=13;  % which day of the electricity file

P.m=10000;              % number of vehicles
P.modechoice=false;     % mode choice?
P.carbonprice=0;        % $/ton

% simulation options
Sim.e=2;                  % time step length in minutes
Sim.mpcpredict=true;      % perfect prediction?

% technical parameters
Tech.battery=50;                % kWh
Tech.chargekw=20;               % kW
Tech.consumption=0.15*30/60;    % kWh/min   (0.15 kwh/km * 30 km/h / 60 min/h)
Tech.cyclingcost=100/4000;%200/2000;      % $/kWh     (batterycost[$/kWh] / lifelength[cycles])
Tech.efficiency=0.9;

% operational parameters
Operations.initialsoc=0.7;
Operations.minsoc=0.2;
Operations.maxsoc=1;
Operations.v2g=true;
Operations.v2gminsoc=0.5;
Operations.maxwait=Inf; % minutes

% % transport layer: opti
% Relocation.alg='opti'
% Relocation.thor=8;           % horizon (time steps)
% Relocation.rho1=0.01;        % weight of secondary objective
% Relocation.rho2=0.01;        % weight of charging objective for electricity price
% Relocation.rho3=0.01;        % weight of charging objective for SOC
% Relocation.rho4=0.000001;    % weight for fixed charge

% transport layer: simplified
Relocation.alg='simplified';
Relocation.tx=10;  % minutes
Relocation.ts=24; 
Relocation.tr=20; 
Relocation.bmin=0;

% energy layer: aggregate
Charging.mthor=1440;      % macro horizon (minutes)
Charging.extrasoc=0.25; %0.25 % extra soc for energy layer to account for aggregate uncertainty
Charging.beta=30;              % frequency of energy layer simulation (in minutes)

% append to main struct
P.Sim=Sim;
P.Tech=Tech;
P.Operations=Operations;
P.Charging=Charging;
P.Relocation=Relocation;

clear Sim Tech Operations Charging Relocation

% % trip file associated with scenario & case-dependent adjustments
% switch Scenario 
%     case 'NYC2016'
%         % P.tripfile='NYC2016_Jan13-Mar16_10days';
%         % P.tripfolder='NYC2016';
%         P.tripfolder='NYC2016-nodes';
%         P.Operations.maxwait=Inf; % minutes
%         P.tripday=13;     % which day of the trip file
%     case 'NYC2016-10clusters'
%         P.Operations.maxidle=5; % minutes
%         P.tripfolder='NYC2016-nodes';
%         P.Operations.maxwait=Inf; % minutes
%         P.tripday=13;     % which day of the trip file
%     case 'NYC2016-20clusters'
%         P.Operations.maxidle=5; % minutes
%         P.tripfolder='NYC2016-nodes';
%         P.Operations.maxwait=Inf; % minutes
%         P.tripday=13;     % which day of the trip file
%     case 'NYC2016-50clusters'
%         P.Operations.maxidle=5; % minutes
%         P.tripfolder='NYC2016-nodes';
%         P.Operations.maxwait=Inf; % minutes
%         P.tripday=13;     % which day of the trip file
%     case 'NYC2018'
%         % P.tripfolder='NYC2018_10wed';
%         % P.tripfile='NY_trips_10wed_0103-0307_minutes';
%         P.tripfolder='NYC2018_10wed-nodes';
%         P.tripday=1;     % which day of the trip file
%     case 'NYC2016-small'
%         P.tripfile='NYC2016-small_13Jan';
% %         P.tripfolder='NYC2016-small_100';
%         P.tripday=1;     % which day of the trip file
%         P.m=30; 
%     case 'NYC2016-small2'
%         P.scenario='NYC2016-small';
%         P.tripfile='NYC2016-small2_13Jan';
%         P.m=5;
%     case 'Tokyo189'
%         P.tripfile='Tokyo2008_1day_48k-nodes';
%         P.m=3500;
%     case 'Munich'
%         P.tripfile='Munich';
%     case 'Munich_clustered'
%         P.tripfile='Munich_clustered';
%         P.gridfile='Germany_DA_2019';
%         P.Operations.maxwait=Inf;
%         P.Tech.chargekw=22;
%         P.m=5000;
%         P.e=1;
%     case 'Munich_clustered_week'
%         P.tripfolder='Munich_1week';
%         P.gridfile='Germany_DA_2019';
%         P.Operations.maxwait=Inf;
%         P.Tech.chargekw=22;
%         P.m=5000;
%         P.e=1;
%         P.tripday=1;
%         P.gridday=1;
% end


