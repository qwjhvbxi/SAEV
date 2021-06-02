%% Z=CHARGINGMODULE(Par,q,travelMinutes,electricityPrice,co2Emissions)
% Finds optimal aggregate fleet charging set points.
% 
% == Input == 
% Par: struct with parameters:
%   Beta        [scalar]   length of a time step in energy layer
%   chargingHorizon [scalar]   number of time steps in energy layer
%   consumption [scalar]   vehicle consumption (kWh/min)
%   chargekw    [scalar]   vehicle max. charging power (kW)
%   battery     [scalar]   vehicle battery capacity (kWh)
%   minsocfleet [scalar]   minimum state of charge in aggregate
%   maxsoc      [scalar]   maximum state of charge
%   efficiency  [scalar]   battery round-trip charging efficiency 
%   carbonprice [scalar]   carbon price ($/ton)
%   cyclingcost [scalar]   equivalent cost of 1 battery cycle
%   v2g         [scalar logical]   use V2G (discharge to the grid)?
% q                 [m x 1]    state of charge of each vehicle
% travelMinutes     [chargingHorizon x 1] aggregate minutes to travel  
% electricityPrice  [chargingHorizon x 1] electricity price ($/kWh) 
% co2Emissions      [chargingHorizon x 1] co2 emissions from grid (g/kWh) 
% 
% == Output ==
% Z: [4 x 1]    charging set points and info until next call, in the form: 
%               [relative charging, relative discharging, max charging allowed, energy required by trips]
% 
% See also: mainsim

function Z=chargingmodule(Par,q,travelMinutes,electricityPrice,co2Emissions)

m=length(q);

% energy layer variable: static values
E.v2g=Par.v2g;                  % use V2G?
E.efficiency=Par.efficiency;    % roundtrip (discharge) efficiency
E.socboost=1e4;                 % soft constraint for final soc (TODO: change depending on inputs)
E.T=Par.chargingHorizon;        % number of time steps in energy layer
E.cyclingcost=Par.cyclingcost;  % battery cycling cost [$/kWh]
E.storagemax=Par.battery*m*Par.maxsoc;    % max total energy in batteries [kWh]
E.carbonprice=Par.carbonprice;                % carbon price [$ per kg]
E.maxchargeminute=Par.chargekw/60*Par.aggregateratio;    % energy exchangeable per minute per vehicle [kWh]

actualminsoc=min(Par.minsocfleet,mean(q)*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of energy layer

E.storagemin=Par.battery*Par.m*actualminsoc; % kWh
E.einit=sum(q)*Par.battery;     % total initial energy [kWh]
E.etrip=travelMinutes*Par.consumption;   % energy used per step [kWh] 

E.dkav=max(0,m*Par.Beta-travelMinutes);         % minutes of availability of cars
E.electricityprice=electricityPrice; % convert to [$/kWh]
E.emissionsGridProfile=co2Emissions; % [g/kWh]

maxc=E.dkav*E.maxchargeminute; % max exchangeable energy per time step [kWh]

ELayerResults=aevopti11(E);

if ~isempty(ELayerResults)

    % Z: [relative charging, relative discharging, max charging allowed, energy required by trips]
    Z=[ ELayerResults.charging(1) , ELayerResults.discharging(1)*logical(Par.v2g) , maxc(1) , E.etrip(1) ]';
    Z(isnan(Z))=0;

else

    % in case there is no feasible solution, charge as much as possible
    Z=[max(1,maxc(1)),0,max(1,maxc(1)),E.etrip(1)]';

end

end


