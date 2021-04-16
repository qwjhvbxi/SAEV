%% [e,ef]=setpointvehicle(Par,q,f,s,SetPoints)
% 

function [e,ef]=setpointvehicle(Par,q,s,f,SetPoints)

% power exchanged for vehicles charging
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

% needed FCR power
FCRNeed=(f-50)/(Par.limits(2)-Par.limits(1))*2;
FCRNeedUp=Par.af*min(1,max(0,FCRNeed)); % charge
FCRNeedDown=Par.af*min(1,max(0,-FCRNeed)); % discharge

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

end