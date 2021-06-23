%% [SetPoints]=SETPOINTFLEET(Par,q,s,z)
% Create set points at beginning of charging period.
% 
% See also: mainsim, chargingsetpoints, setpointvehicle

function [SetPoints]=setpointfleet(Par,q,s,z)

H=Par.Beta/Par.Epsilon;

ac=Par.chargekw/Par.battery/60*Par.Epsilon;    % charge rate per time step (normalized)

% power exchanged for vehicles charging
acv=(q<Par.fastchargesoc)*ac+(q>=Par.fastchargesoc)*ac*Par.slowchargeratio;

% aggregate set point (kWh)
SetPointUpPeriod=(z(1));  % set point of aggregate fleet (kWh)
SetPointDownPeriod=(z(2));  % set point of aggregate fleet (kWh)

% expected total vehicle capacity in the period (kWh)
CapUpPeriod=max(0,s.*min(acv*H,Par.maxsoc-q)*Par.battery); % charge
CapDownPeriod=max(0,s.*min(acv*H,(q-Par.v2gminsoc)*Par.efficiency)*Par.battery); % discharge

% set point for each vehicle for each time step (kWh)
SetPoints(1)=min(SetPointUpPeriod,sum(CapUpPeriod))/H;  % set point of aggregate fleet (kWh) UP
SetPoints(2)=min(SetPointDownPeriod,sum(CapDownPeriod))/H;  % set point of aggregate fleet (kWh) DOWN

end