%% CPAR
% Create default P struct.
% 
% see also SCENARIOGENERATOR


% scenario data
P.scenario='NYC2016';
P.tripfile='NYC2016_Jan13-Mar16_10days';
P.scenarioid=1;
% load(['data/scenarios/' P.scenario{1} '.mat'],'T');
load(['data/scenarios/' P.scenario '.mat'],'T','C');
P.coords=C;
P.T=T;

% electricity profile 
P.eleproftype='TokyoDA-FY2017-Reduced';
P.eleprofseed=311;%u(1);
load(['data/' P.eleproftype '.mat'],'u','x');
P.melep=repmat(x(:,P.eleprofseed),2,1);      % macro elep

% clear T u x

P.chargingstations=[4,4;5,6];



% model options
P.mpcpredict=1;
P.enlayeralg='aggregate';%'opti';
P.trlayeralg='simplified';%'opti';
P.v2g=true;

% model parameters: times
P.tau=1;            % dataset resolution scaling factor
P.e=2;              % time step length in minutes
P.beta=15;
P.thor=10;          % horizon (time steps)
P.mthor=48;         % macro horizon (macro time step)

% model parameters: weights
P.rho1=0.01;        % weight of secondary objective
P.rho2=0.01;        % weight of charging objective for electricity price
P.rho3=0.01;        % weight of charging objective for SOC
P.rho4=0.000001;    % weight for fixed charge

% model parameters: fleet
P.m=10000;             % number of vehicles
P.battery=50;
P.chargekw=20;
P.consumption=0.15*30/60;   % consumption per minute: 0.15 kwh/km * 30 km/h / 60 min/h
P.cyclingcost=20000/2000;   % batterycost[yen/kWh] / lifelength[cycles] -> 10 yen/kWh

% operational parameters
P.initialsoc=0.7;
P.minsoc=0.2;
P.maxsoc=1;
P.mminsoc=0.3;  % macro min soc
P.mmaxsoc=1;    % macro max soc
P.v2gminsoc=0.5;
P.maxwait=10;

% model parameters: aggregate
P.tx=5;
P.ts=15;
P.tr=15;
P.bmin=0;