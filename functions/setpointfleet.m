%% [SetPoints]=SETPOINTFLEET(Par,q,s,u,z)
% Create set points (kWh) at beginning of charging period.
% 
% See also: mainsim, chargingsetpoints, setpointvehicle

function [SetPoints]=setpointfleet(Par,q,s,u,z)

% number of time steps
H=Par.SPlength/Par.Epsilon;

% number of charging stations
ncs=length(Par.cssize);

% charge rate per time step (normalized)
ac=Par.chargekw/Par.battery/60*Par.Epsilon;

% power exchanged for vehicles charging
acv=(q<Par.fastchargesoc)*ac+(q>=Par.fastchargesoc)*ac*Par.slowchargeratio;

% aggregate set point (kWh)
SetPointUpPeriod=(z(1)/Par.Beta*Par.SPlength);  % set point of aggregate fleet (kWh)
SetPointDownPeriod=(z(2)/Par.Beta*Par.SPlength);  % set point of aggregate fleet (kWh)

% expected capacity in the period for each vehicle (kWh)
CapUpPeriod=max(0,s.*min(acv*H,Par.maxsoc-q)*Par.battery); % charge
CapDownPeriod=max(0,s.*min(acv*H,(q-Par.v2gminsoc)*Par.efficiency)*Par.battery); % discharge

% expected capacity in the period for each charging station (kWh)
CapUpPeriodCS=min(accumarray(u(s)',CapUpPeriod(s)',[ncs,1]),Par.cssize*Par.chargekw*Par.SPlength/60);
CapDownPeriodCS=min(accumarray(u(s)',CapDownPeriod(s)',[ncs,1]),Par.cssize*Par.chargekw*Par.SPlength/60);

% set point for fleet for each time step (kWh)
SetPoints(1)=min(SetPointUpPeriod,sum(CapUpPeriodCS))/H;  % set point of aggregate fleet (kWh) UP
SetPoints(2)=min(SetPointDownPeriod,sum(CapDownPeriodCS))/H;  % set point of aggregate fleet (kWh) DOWN

end



