%% [e,ef]=SETPOINTVEHICLE(Par,q,s,u,SetPoints,f)
% Return set point for each vehicle given fleet set points 'SetPoints'
% 
% See also: mainsim, setpointfleet

function [e,ef]=setpointvehicle(Par,q,s,u,SetPoints,f)

if nargin<6
    f=0;
end

% charge rate per time step (normalized)
Par.ac=Par.chargekw/Par.battery/60*Par.Epsilon;    


%% select which vehicles are actually connected

% number of charging stations (CS)
ncs=length(Par.cssize);

% number of vehicles at each CS
Occupancy=histc(u,1:ncs);

% find overcrowded CS
Overcrowded=find(Occupancy>Par.cssize');
for i=1:length(Overcrowded)
    
    % find vehicles at overcrowded CS
    thisStation=Overcrowded(i);
    atThisStation=find(u==thisStation);
    
    % only soc extremes are connected. Set capacity to zero for vehicles that are not available
    [~,I]=sort(q(atThisStation));

    if SetPoints(2)>0
        % for discharging, select highest soc vehicles at station within the limit
        s(atThisStation(I(1:end-Par.cssize(thisStation))))=0;
    else
        % for charging, select lowest soc vehicles at station within the limit
        s(atThisStation(I(Par.cssize(thisStation)+1:end)))=0;
    end
    
end


%% charging

% max power exchangeable for vehicles charging
acv=(q<Par.fastchargesoc)*Par.ac+(q>=Par.fastchargesoc)*Par.ac*Par.slowchargeratio;

% available power from fleet
maxsoceff=1; % maxsoceff=Par.maxsoc;
v2gallowed=q>Par.v2gminsoc;

% actual capacity in this time step
CapUp=s.*min(acv,maxsoceff-q); % charge
CapDown=s.*v2gallowed.*min(acv,(q-Par.minsoc)*Par.efficiency); % discharge

% calculate ratios
if sum(CapUp)>0
    eRatioUp=min(1,SetPoints(1)/(sum(CapUp)*Par.battery));
else
    eRatioUp=0;
end
if sum(CapDown)>0
    eRatioDown=min(1,SetPoints(2)/(sum(CapDown)*Par.battery));
else
    eRatioDown=0;
end

% calculate charging for each vehicle
e=CapUp*eRatioUp-CapDown*eRatioDown;


%% FCR provision

if f~=0 && Par.fcrcontracted>0

    af=Par.fcrcontracted*1000/Par.battery/60*Par.Epsilon;    % FCR rate per time step (normalized)

    % needed FCR power
    FCRNeed=(f-50)/(Par.fcrlimits(2)-Par.fcrlimits(1))*2;
    FCRNeedUp=af*min(1,max(0,FCRNeed)); % charge
    FCRNeedDown=af*min(1,max(0,-FCRNeed)); % discharge

    % available power from fleet
    AvailableUp=s.*min(Par.ac-max(0,e),1-(q+e)); % charge
    AvailableDown=s.*min(Par.ac-max(0,-e),(q+e)-Par.minsoc); % discharge

    % calculate ratios
    if sum(AvailableUp)>0
        FCRRatioUp=min(1,FCRNeedUp/sum(AvailableUp));
    else
        FCRRatioUp=0;
    end
    if sum(AvailableDown)>0
        FCRRatioDown=min(1,FCRNeedDown/sum(AvailableDown));
    else
        FCRRatioDown=0;
    end

    % calculate FCR power contributions
    ef=AvailableUp.*FCRRatioUp-AvailableDown.*FCRRatioDown;

else
    
    ef=0;
    
end

end