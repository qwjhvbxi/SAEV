% Vin: vehicles information in the form: [station delay soc charging relocating]


function [V,b]=relocationmodule(Vin,Par)

nc=size(Par.Trs);
nv=size(Vin,1);

position=Vin(:,1);
delay=Vin(:,2);
soc=Vin(:,3);
charging=Vin(:,4);
relocating=Vin(:,5);
used=zeros(nv,1);
if isfield(Par,'LimitFCR')
    LimitFCR=Par.LimitFCR;
else
    LimitFCR=0;
end

% vehicles at clusters
uci=position.*(delay==0); % idle
ucr=position.*relocating; % relocating

% number of vehicles at each station
uv=histc(uci,1:nc);

% vehicles relocating here between now and now+ts 
uvr=histc(ucr,1:nc);

% expected imbalance at stations
b=uv-Par.dw+Par.a_ts-Par.a_to+uvr;

% identify feeder and receiver stations
F=min(uv,(b-Par.bmin).*(b>=Par.bmin)); % feeders
R=(-b+Par.bmin).*(b<Par.bmin); % receivers

% identify optimal relocation flux
x=optimalrelocationfluxes(F',R',Par.Trs,Par.limite);

if ~isempty(x)

    % read results
    [Fs,Rs,Vr]=find(x);

    % distance of relocation
    ReloDistance=Par.Trs(sub2ind(size(Par.Trs),Fs,Rs)); 

    % duplicate fluxes with multiple vehicles
    Fs=repelem(Fs,Vr);
    Rs=repelem(Rs,Vr);
    ReloDistance=repelem(ReloDistance,Vr);

    % satisfy longer relocation tasks first
    [ReloDistanceSorted,dstnid]=sort(ReloDistance,'descend');

    % for each single task (1 vehicle)
    for ka=1:length(ReloDistanceSorted)

        % find candidate vehicles for the task with enough soc and idle, give
        % priority to non-charging vehicles
        candidates=(1+soc-charging*Par.chargepenalty).*(position==Fs(dstnid(ka))).* ...
            (soc/Par.ad >= delay+ReloDistanceSorted(ka)).*(relocating==0).*(delay+ReloDistanceSorted(ka)<=Par.limite); 

        % remove unavailable vehicles
        candidates(candidates==0)=NaN;
        
        % remove charging vehicles if exceeding the LimitFCR
        if sum(candidates<0)<=LimitFCR
            candidates(candidates<0)=NaN;
        end

        % sort candidate vehicles by SOC
        [ur,ui]=max(candidates);

        % if there is a vehicle available
        if ~isnan(ur)

            % update destination station
            position(ui)=Rs(dstnid(ka));

            % update status
            used(ui)=true;
            relocating(ui)=1;

        end
    end
    
    V=[position , used];
    
else
    
    V=[Vin(:,1) , zeros(nv,1)];
    
end

end