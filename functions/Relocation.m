% Vin: vehicles information in the form: [station delay soc charging relocating]


function [V,b]=Relocation(Vin,Par)

nc=size(Par.Trs);
nv=size(Vin,1);

position=Vin(:,1);
delay=Vin(:,2);
soc=Vin(:,3);
relocating=Vin(:,5);
used=zeros(nv,1);

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
x=optimalrelocationfluxes(F',R',Par.Trs);

if ~isempty(x)

    % read results
    [Fs,Rs,Vr]=find(x);

    % distance of relocation
    arri=Par.Trs(sub2ind(size(Par.Trs),Fs,Rs)); 

    % duplicate fluxes with multiple vehicles
    Fs=repelem(Fs,Vr);
    Rs=repelem(Rs,Vr);
    arri=repelem(arri,Vr);

    % satisfy longer relocation tasks first
    [arris,dstnid]=sort(arri,'descend');

    for ka=1:length(arris)

        % find candidate vehicles for the task with enough soc and idle
        candidates=soc.*(position==Fs(dstnid(ka))).*(soc/Par.ad >= arris(ka)).*(relocating==0); 

        % remove unavailable vehicles
        candidates(candidates==0)=NaN;

        % sort candidate vehicles by SOC
        [ur,ui]=max(candidates);

        % if there is a vehicle available
        if ~isnan(ur)

            % update destination station
            position(ui)=Rs(dstnid(ka));

            % update status
            used(ui)=true;

        end
    end
    
    V=[position , used];
    
else
    
    V=[Vin(:,1) , zeros(nv,1)];
    
end

end