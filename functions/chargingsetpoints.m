%% [e,ef]=CHARGINGSETPOINTS(Par,qi,s,Z[,setPoints,f])
% Finds charging set points for each vehicle.
% 
% == Input == 
% Par: struct with parameters:
%   Beta        [scalar]   length of a time step in charging module (minutes)
%   Epsilon     [scalar]   length of a time step in simulation (minutes)
%   consumption [scalar]   vehicle consumption (kWh/minute)
%   chargekw    [scalar]   vehicle max. charging power (kW)
%   battery     [scalar]   vehicle battery capacity (kWh)
%   minsoc      [scalar]   minimum state of charge 
%   maxsoc      [scalar]   maximum state of charge
%   efficiency  [scalar]   battery round-trip charging efficiency 
%   v2gminsoc   [scalar]   minimum state of charge to participate in V2G discharging
%   csp         [scalar logical]   use constant setpoint? 
%  --- if csp==0: ---
%   [optional]  refillmaxsoc    [scalar]  state of charge below which vehicles will charge in any case
%  --- if csp==1: ---
%   [optional]  fastchargesoc   [scalar]  state of charge below which vehicles will charge at full power
%   [optional]  slowchargeratio [scalar]  ratio of full power at which vehicles will charge when soc is abobe fastchargesoc
%   [optional]  fcrcontracted   [scalar]  FCR contracted power (MW)
%   [optional]  fcrlimits       [scalar]  FCR frequency limits
% q: [m x 1]    state of charge of each vehicle
% s: [m x 1]    status of each vehicle
% Z: [4 x 1]    charging set points and info until next call, in the form: 
%               [charging, discharging, max charging allowed, energy required by trips]
% [optional]  newSetPoint:  [scalar logical] start new set point? (default: 0)
% [optional]  f:    [scalar]   current grid frequency level (Hz)
% 
% == Output ==
% e:  [m x 1]   charging power for each vehicle as a fraction of SOC
% ef: [m x 1]   charging power for each vehicle for frequency control as a fraction of SOC (when applicable)
% 
% See also: chargingmodule, mainsim

function [e,ef]=chargingsetpoints(Par,q,s,Z,setPoints,f)

if nargin<6
    f=0;
end
if ~isfield(Par,'refillmaxsoc')
    Par.refillmaxsoc=0;
end

Par.ac=Par.chargekw/Par.battery/60*Par.Epsilon;    % charge rate per time step (normalized)
% Par.ad=Par.consumption/Par.battery*Par.Epsilon;    % discharge rate per time step (normalized)

if Par.csp  % setpoint based

%     % default values
%     if ~isfield(Par,'fastchargesoc')
%         Par.fastchargesoc=1;
%     end
%     if ~isfield(Par,'slowchargeratio')
%         Par.slowchargeratio=1;
%     end
%     if ~isfield(Par,'fcrlimits')
%         Par.fcrlimits=[0,0];
%     end
%     if ~isfield(Par,'fcrcontracted')
%         Par.fcrcontracted=0;
%     end    
    
    [e,ef]=setpointvehicle(Par,q,s,f,setPoints);


else     % capacity based
    
    m=length(q);

    v2gallowed=q>Par.v2gminsoc;
    extracharge=(q<Par.refillmaxsoc);
    chargevector=max(-1,min(1,(ones(1,m)*(Z(1)/Z(3))-v2gallowed*(Z(2)/Z(3))+extracharge)))*Par.ac;

    capUp=s.*min(Par.ac,Par.maxsoc-q); % charge
    capDown=s.*min(Par.ac,(q-Par.minsoc)*Par.efficiency); % discharge

    e=min(capUp,max(0,chargevector))+max(-capDown,min(0,chargevector));
    ef=0;

end

end
