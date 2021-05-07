

function zmacro=chargingmodule(Par,qi,dktripnow,melepnow,mco2now)

% energy layer variable: static values
E.v2g=Par.v2g;             % use V2G?
E.efficiency=Par.efficiency;     % roundtrip (discharge) efficiency
E.socboost=1e4;                     % soft constraint for final soc (TODO: change depending on inputs)
E.T=Par.chargingHorizon;                % number of time steps in energy layer
E.cyclingcost=Par.cyclingcost;                       % battery cycling cost [$/kWh]
E.storagemax=Par.battery*Par.m*Par.maxsoc;    % max total energy in batteries [kWh]
E.carbonprice=Par.carbonprice;                            % carbon price [$ per kg]
E.maxchargeminute=Par.ac*Par.battery/Par.e*Par.aggregateratio;    % energy exchangeable per minute per vehicle [kWh]

actualminsoc=min(Par.minsoc+Par.extrasoc,mean(qi)*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of energy layer

E.storagemin=Par.battery*Par.m*actualminsoc; % kWh
E.einit=sum(qi)*Par.battery;     % total initial energy [kWh]
E.etrip=dktripnow*Par.ad*Par.battery/Par.e;   % energy used per step [kWh] 

E.dkav=max(0,Par.m*Par.Beta-dktripnow);         % minutes of availability of cars
E.electricityprice=melepnow; % convert to [$/kWh]
E.emissionsGridProfile=mco2now; % [g/kWh]

maxc=E.dkav*E.maxchargeminute; % max exchangeable energy per time step [kWh]

ELayerResults=aevopti11(E);

if ~isempty(ELayerResults)

    % TODO: should be only first value
    % zmacro: [relative charging, relative discharging, max charging allowed, energy required by trips]
    zmacro=[ ELayerResults.charging(1) , ELayerResults.discharging(1)*logical(Par.v2g) , maxc(1) , E.etrip(1) ]';
    zmacro(isnan(zmacro))=0;

else

    % in case there is no feasible solution, charge as much as possible
    zmacro=[1,0,maxc(1),E.etrip(1)]';

end

end


