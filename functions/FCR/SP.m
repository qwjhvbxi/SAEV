%% [SetPoints]=SP(Par,q,s,z)
% create set points at beginning of charging period

function [SetPoints]=SP(Par,q,s,z)

% power exchanged for vehicles charging
acv=(q<Par.fastchargesoc)*Par.ac+(q>=Par.fastchargesoc)*Par.ac*Par.slowchargeratio;

% aggregate set point (kWh)
SetPointUpPeriod=(z(1));  % set point of aggregate fleet (kWh)
SetPointDownPeriod=(z(2));  % set point of aggregate fleet (kWh)

% expected total vehicle capacity in the period (kWh)
CapUpPeriod=max(0,s.*min(acv*Par.H,Par.maxsoc-q)*Par.battery); % charge
CapDownPeriod=max(0,s.*min(acv*Par.H,(q-Par.v2gminsoc)*Par.efficiency)*Par.battery); % discharge

% set point for each vehicle for each time step (kWh)
SetPoints(1)=min(SetPointUpPeriod,sum(CapUpPeriod))/Par.H;  % set point of aggregate fleet (kWh) UP
SetPoints(2)=min(SetPointDownPeriod,sum(CapDownPeriod))/Par.H;  % set point of aggregate fleet (kWh) DOWN

end