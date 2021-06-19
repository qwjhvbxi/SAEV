%% [Res]=MAINSIM(P,dispiter)
% Main simulation.
% Vehicles start at beginning of time step, arrive at end of time step.
% 
% Problem: pricing cannot be run with charging because of how EMD is
% calculated right now
% Problem: OD-based pricing cannot run with aggregate prediction
% 
% See also: main

function [Res]=mainsim(P,dispiter)

DataFolder=getdatafolder();

%% load external files: scenario, trips 

% load distance matrix
[T,D,clusters,clusterIDs,chargingStations]=getscenario(P.scenario);

% load trips
[A,Atimes,cumulativeTripArrivals,~]=gettrips(P);

% load electricity prices and carbon emissions
[elepMinute,co2Minute,~]=readexternalfile([DataFolder 'grid/' P.gridfile '.csv'],P.gridday,true);


%% setup clustering

n=size(clusters,1);             % number of nodes 
As=clusters(A(:,1:2));          % OD at clusters
nc=length(clusterIDs);    % number of clusters


%% load predictions

% Aforecast is a t x n^2 matrix with probabilities of each OD pair at each
% time interval (needed only for pricing optimization)
[fo,fd,dkod,Aforecast]=loadpredictions(P,As,Atimes);


%% initialize parameters of simulation

% parameters
r=cumulativeTripArrivals(1441);           % number of requests
tsim=1440/P.Sim.e;        % number of time steps

% electricity and emissions profiles
elep=average2(elepMinute,P.Sim.e);
co2=average2(co2Minute,P.Sim.e);
f=sparse(tsim,1);

% main simulation variables
q=zeros(tsim,P.m);          % SOC
u=zeros(tsim,P.m,'uint16'); % vehicles in charging stations
d=zeros(tsim,P.m,'uint16'); % delay
e=zeros(tsim,P.m);          % charging
ef=zeros(tsim,P.m);         % FCR charging

% working variables
g=zeros(1,P.m);             % current idle time
s=false(5,P.m);             % current status: [relocating, connected, moving to charging station]
queue=zeros(100,1);         % temporary variable to store queued arrivals

% results variables
waitingestimated=zeros(r,1);    % estimated minutes to wait for each request
waiting=zeros(r,1);             % minutes waited for each request
dropped=false(r,1);             % request is dropped?
chosenmode=false(r,1);          % which mode is chosen?
offeredprices=zeros(r,1);       % price offered to each passenger
status=zeros(tsim,P.m);         % vehicle status
relodist=zeros(tsim,1);         % distances of relocation (at moment of decision)
relodistkm=zeros(tsim,1);       % distances of relocation in km (at moment of decision)
tripdist=zeros(tsim,1);         % distances of trips (at moment of acceptance)
tripdistkm=zeros(tsim,1);       % distances of trips in km (at moment of acceptance)
% pooling=zeros(r,1);             % pool ID of each user (if ride shared)


%% setup internal parameters

Par=struct('D',D,'Epsilon',P.Sim.e,'minsoc',P.Operations.minsoc,'maxsoc',P.Operations.maxsoc,'modechoice',P.modechoice,...
    'battery',P.Tech.battery,'maxwait',P.Operations.maxwait,'VOT',P.Pricing.VOT,'WaitingCostToggle',P.Pricing.pricingwaiting,...
    'LimitFCR',0,'chargepenalty',1,'v2gminsoc',P.Operations.v2gminsoc,'efficiency',P.Tech.efficiency,...
    'csp',false,'refillmaxsoc',0,'aggregateratio',1,'chargekw',P.Tech.chargekw,'consumption',P.Tech.consumption);

% distances
if ~isstruct(T)
    Tr=max(1,round(T/P.Sim.e));% distance matrix in steps
    Tr(1:n+1:end)=0;          % no distance between same node
    Trs=Tr(clusterIDs,clusterIDs);
    Par.Tr=Tr;
    [~,closestCS]=min(Tr(:,chargingStations(:,1)),[],2); % closest charging station to each node
end
tripDistances=nan(r,1);
tripDistancesKm=D(sub2ind(size(D),A(:,1),A(:,2)))/1000;


%% setup relocation module

if isfield(P,'Relocation') && ~isempty(P.Relocation)
    ts=round(P.Relocation.ts/P.Sim.e);
    tr=round(P.Relocation.tr/P.Sim.e);
    tx=round(P.Relocation.tx/P.Sim.e);
    bmin=P.Relocation.bmin;
    b=zeros(ceil(tsim/tx),nc,'double');             % imbalance
    autoRelocation=1;
else
    b=NaN;
    autoRelocation=0;
end


%% setup frequency control reserve module

if isfield(P,'Charging') && isfield(P,'FCR') && ~isempty(P.FCR) 
    
    Par.csp=true;
    if isfield(P.FCR,'aggregatechargeratio')
        Par.aggregateratio=P.FCR.aggregatechargeratio; % charge rate at aggregate level optimization
    end
    
    [fraw,~,fresolution]=readexternalfile([DataFolder 'grid/' P.FCR.filename],P.gridday,false);
    f=average2(fraw,P.Sim.e/fresolution);
    
    Par.LimitFCR=ceil(P.FCR.contracted*1000/P.Tech.chargekw);
    Par.fcrcontracted=P.FCR.contracted;
    Par.fcrlimits=P.FCR.limits;
    Par.fastchargesoc=P.FCR.fastchargesoc;
    Par.slowchargeratio=P.FCR.slowchargeratio;
    
end


%% setup charging module

dynamicCharging=false;
Trips=[];
Beta=0;

if isfield(P,'Charging') && ~isempty(P.Charging)
    
    if strcmp(P.Charging,'night')
        
        % night charging
        limitHour=5;        % charging hour limit
        Par.refillmaxsoc=0.6;   % min. soc to trigger charging during day
        zmacro=[ [1;0;1;0]*ones(1,60*limitHour/P.Sim.e) , [0;0;1;0]*ones(1,60*(24-limitHour)/P.Sim.e) ];
        
    else
    
        dynamicCharging=true;
        
        % inputs
        Beta=P.Charging.beta;
        chargingHorizon=round(P.Charging.mthor/Beta);
        
        % variables
        etsim=floor(1440/Beta); % number of charging decisions
        melep=average2(elepMinute,Beta);
        mco2=average2(co2Minute,Beta);

        % generate aggregate trip statistics
        emdFileName=[DataFolder 'temp/emd-' P.tripfolder '-' num2str(P.tripday) '-' num2str(Beta) '.mat'];
        if exist(emdFileName,'file')
            load(emdFileName,'dkemd','dkod','dktrip');
        else
            dkod=computetraveltime(A,Atimes,T,Beta);% TODO: should be on prediction
            dkemd=computeemd(fo,fd,T,Beta,clusterIDs); 
            dktrip=dkod+dkemd;
            save(emdFileName,'dkemd','dkod','dktrip');
        end
        Trips=struct('dkod',dkod,'dkemd',dkemd,'dktrip',dktrip);
        
        % energy layer variable: static values
        Par.Beta=P.Charging.beta;
        Par.chargingHorizon=round(P.Charging.mthor/Par.Beta);
        Par.v2g=P.Operations.v2g;
        Par.cyclingcost=P.Tech.cyclingcost;
        Par.m=P.m;
        Par.carbonprice=P.carbonprice;
        Par.minsocfleet=Par.minsoc+P.Charging.extrasoc;
        
        % matrix of optimal control variables for energy layer
        zmacro=zeros(4,etsim+chargingHorizon); 
        
    end

else
    
    % charge as much as possible
    zmacro=[1;0;1;0]*ones(1,tsim);
    
end


%% setup pricing module 

addpath functions/pricing

% add info to Pricing struct
P.Pricing.relocation=autoRelocation;
P.Pricing.c=D(clusterIDs,clusterIDs)/1000;

% initializations
perDistanceTariff=ones(nc,nc).*P.Pricing.basetariffkm; % matrix of fare per minute
surchargeMat=zeros(nc,nc);  % surcharges per stations

if isempty(P.Pricing.alternativecost)
    % alternative price for each user calculated with trip distances in km
    Aaltp=P.Pricing.alternativecostkm*tripDistancesKm; 
else
    Aaltp=P.Pricing.alternativecost;
end

altp=0;

if P.Pricing.dynamic
    tp=round(P.Pricing.tp/P.Sim.e);       % pricing interval
    tpH=round(P.Pricing.horizon/P.Sim.e); % pricing horizon
else
    tpH=0;
    tp=tsim;
end

% initialization for dynamic pricing
tariff=ones(nc^2,ceil(tsim/tp))*P.Pricing.basetariffkm;
surcharge=zeros(nc*2,ceil(tsim/tp));


%% initial states

q(1,:)=P.Operations.initialsoc.*ones(1,P.m);      % initial state of charge
if isfield(P.Operations,'uinit')
    u(1,:)=P.Operations.uinit;
else
    u(1,:)=chargingStations(randi(nc,1,P.m),1);                 % initial position of vehicles
end
if isfield(P.Operations,'dinit')
    d(1,:)=P.Operations.dinit;
end

% initial status
atChargingStation=sum(u(1,:)==chargingStations(:,1));
s(1,:)=logical(atChargingStation.*(d(1,:)==0));
s(2,:)=logical(~atChargingStation.*(d(1,:)==0));

% setup for trip generation
cumulativeTripArrivals=cumulativeTripArrivals(1:P.Sim.e:end);


%% variables for progress display and display initializations

S.starttime=cputime;
S.lasttime=S.starttime;
S.comptime=[cputime;zeros(tsim,1)];


%% start of iterations

for i=1:tsim
    
    
	%% display progress
    
    displayprogress(i,tsim,dispiter,S.comptime(i)-S.comptime(1),40)
    ui=double(u(i,:));
    di=double(d(i,:));
    relodist(i)=0;
    relodistkm(i)=0;
    
    
    %% calculate current distance matrix
    
    if isstruct(T)
        [thisT]=gettraveltimenow(T,i*Par.Epsilon);
        Tr=max(1,round(thisT/P.Sim.e));% distance matrix in steps
        Trs=Tr(clusterIDs,clusterIDs);
        Par.Tr=Tr;
        [~,closestCS]=min(Tr(:,chargingStations(:,1)),[],2); % closest charging station to each node
    end
    
    
    %% move idle vehicles back to charging stations
    
    if nc<n
        idleTime=g.*(1-atChargingStation);
        IdleReached=(idleTime>=P.Operations.maxidle/P.Sim.e);
        ui(IdleReached)=chargingStations(closestCS(ui(IdleReached)),1);
        relodistCS=Tr(sub2ind(size(Tr),u(i,IdleReached),ui(IdleReached)));
        relodistCSkm=D(sub2ind(size(D),u(i,IdleReached),ui(IdleReached)))/1000;
        di(IdleReached)=relodistCS;
        s(2,IdleReached)=0;
        s(5,IdleReached)=1;
        relodist(i)=relodist(i)+sum(relodistCS);
        relodistkm(i)=relodistkm(i)+sum(relodistCSkm);
    end
    
    
    %% pricing optimization

    if P.modechoice
        
        if ~P.Pricing.dynamic || mod(i-1,tp)==0
            
            % current pricing number
            kp=ceil(i/tp);
            
            % expected trips
            selection0=cumulativeTripArrivals(i)+1:cumulativeTripArrivals(min(length(cumulativeTripArrivals),i+tpH+1));

            % TODO / PRICING: change pricing with per km? Can be calculated from the start and
            %       they do not change
            % TODO / PRICING: reorganize for imperfect prediction!
            %       price of alternative option should be dependent on OD, not single
            %       passenger (the one seen by optimization). Specific cost only for
            %       mode choice module
            
            AsNow=As(selection0,:);
            altpNow=Aaltp(selection0);

            [perDistanceTariff,surchargeNodes,altp]=pricingmodule(P.Pricing,AsNow,altpNow,clusters(ui));

            tariff(:,kp)=perDistanceTariff(:);
            surcharge(:,kp)=surchargeNodes;

            surchargeMat=surchargeNodes(1:nc)+surchargeNodes(nc+1:2*nc)';
            
        end

        if autoRelocation
        
            option1=exp(-perDistanceTariff.*P.Pricing.c-surchargeMat);
            multiplier=option1./(option1+exp(-altp));

            % expected OD matrices for different future horizons
            % TODO / PRICING: remove reshape, treat all pricing variables as reshaped (from
            % pricingmodule)
            a_ts=multiplier.*reshape(sum(Aforecast(i:i+ts,:)),nc,nc);
            a_to=multiplier.*reshape(sum(Aforecast(i:i+ts+tr,:)),nc,nc);

            % update expected arrivals and departures
            fd(i:i+ts,:)=ones(ts+1,1)*sum(a_ts)/(ts+1);     % arrivals 
            fo(i:i+ts+tr,:)=ones(ts+tr+1,1)*sum(a_to,2)'/(ts+tr+1); % departures
            
        end
        
    end
    

    %% relocation

    % if it's time for a relocation decision
    if autoRelocation && mod(i-1,tx)==0
        
        % current relocation number
        kt=(i-1)/tx+1;

        % number of waiting passenger at station
        dw=histcounts(As(queue(queue>0),1),1:nc+1)';
        
        available=sum(s(1:2,:))';
        
        % Vin: vehicles information in the form: [station delay soc connected relocating]
        Vin=[clusters(ui) , Par.Epsilon*di' , available.*q(i,:)' , s(1,:)' , logical(s(4,:)+s(5,:))' ];
        ParRel.dw=dw; % number of passengers waiting at each station
        ParRel.a_ts=round(sum(fd(i:i+ts,:)))'; % expected arrivals between now and now+ts
        ParRel.a_to=round(sum(fo(i:i+ts+tr,:)))'; % expected requests between now and now+ts+tr
        ParRel.Trs=Trs*Par.Epsilon;
        ParRel.limite=P.Relocation.ts;
        ParRel.bmin=bmin;
        ParRel.consumption=P.Tech.consumption;
        ParRel.battery=P.Tech.battery;
        ParRel.LimitFCR=Par.LimitFCR;
        ParRel.chargepenalty=Par.chargepenalty;
        
        [Vout,bkt]=relocationmodule(Vin,ParRel);
        
        % update vehicles position
        used=logical(Vout(:,2));
        relodestinations=clusterIDs(Vout(used,1));
        relodisti=Tr(sub2ind(size(Tr),ui(used),relodestinations'));
        relodistkmi=D(sub2ind(size(D),ui(used),relodestinations'))/1000;
        
        ui(used)=relodestinations; 
        di(used)=di(used)+relodisti;
        
        % update vehicles status (relocating vehicles cannot be relocated)
        s(1:2,used)=0;
        s(4,used)=1;
        s(5,used)=0;
        
        % update results
        b(kt,:)=bkt;
        relodist(i)=relodist(i)+sum(relodisti);
        relodistkm(i)=relodistkm(i)+sum(relodistkmi);
        
    end


    %% trip assignment

    % generate trip requests for this time step
    trips=(cumulativeTripArrivals(i)+1:cumulativeTripArrivals(i+1))';

    % add previously queued requests and reset queue
    trips=[queue(queue>0);trips];
    queue(:)=0;
    
    if ~isempty(trips)
        
        % calculate pricing    
        selectorClusters=sub2ind(size(Trs),As(trips,1),As(trips,2));
        tripDistances(trips)=Tr(sub2ind(size(Tr),A(trips,1),A(trips,2)))*Par.Epsilon;
        pp=perDistanceTariff(selectorClusters).*tripDistancesKm(trips); 
%             surchargeMat(selectorClusters); % TODO: fix!
        alte=exp(-Aaltp(trips));
        
        % offered prices
        offeredprices(trips)=pp;
        
        % available vehicles depending on status
        available=sum(s(1:3,:))';
        
        % Vin: vehicles information in the form: [station delay soc connected]
        Vin=[ui' , di' , available.*q(i,:)' , s(1,:)' ];
        
        % Bin: passengers info in the form: [O D waiting offeredprice utilityalternative]
        Bin=[A(trips,1:2) , waiting(trips) , pp , alte ];

        if autoRelocation
            [Vout,Bout,tripdisti,tripdistkmi,relodistiPU,relodistkmiPU,queuei]=tripassignmentsaev(Vin,Bin,Par);
        else
            [Vout,Bout,tripdisti,tripdistkmi,relodistiPU,queuei]=tripassignmentcarsharing(Vin,Bin,Par);
            relodistkmiPU=0;
        end
        
        % update vehicles positions
        ui=Vout(:,1)';
        di=Vout(:,2)';
        
        % update vehicles status
        used=logical(Vout(:,3));
        s(1:2,used)=0;
        s(3,used)=1;
        
        % update results
        tripdist(i)=tripdisti;
        tripdistkm(i)=tripdistkmi;
        relodist(i)=relodist(i)+relodistiPU;
        relodistkm(i)=relodistkm(i)+relodistkmiPU;
        chosenmode(trips)=Bout(:,1);
        waiting(trips)=Bout(:,2);
        dropped(trips)=Bout(:,3);
        waitingestimated(trips)=waitingestimated(trips)+Bout(:,4);
        
        queue=trips(queuei(queuei>0));
        
    end
    
    
    %% charging optimization
    
    if dynamicCharging 
        
        if rem(i,Beta/P.Sim.e)==1

            % index of energy layer
            t=(i-1)/(Beta/P.Sim.e)+1;
            
            dktripnow=dktrip(t:t+chargingHorizon-1);    % minutes spent traveling during this horizon
            melepnow=melep(t:t+chargingHorizon-1)/1000; % convert to [$/kWh]
            mco2now=mco2(t:t+chargingHorizon-1); % [g/kWh]

            currentsp=chargingmodule(Par,q(i,:),dktripnow,melepnow,mco2now);
            
            zmacro(:,t)=currentsp;

        end
        
    else
        
        t=i;
        
    end
    
    
    %% simulation variables update
    
    [ei,efi]=chargingsetpoints(Par,q(i,:),s(1,:),zmacro(1:3,t),(rem(i,Beta/P.Sim.e)==1),f(i));
    e(i,:)=ei;
    ef(i,:)=efi;
    
    % update SOC 
    q(i+1,:)=q(i,:)+max(0,e(i,:)+ef(i,:))+min(0,e(i,:)+ef(i,:))/P.Tech.efficiency-(di>0).*(P.Tech.consumption/P.Tech.battery*P.Sim.e);

    % update position
    u(i+1,:)=ui;
    
    % update idle time
    g=(di==0).*(g+1);
    
    % update delay
    d(i+1,:)=max(0,di-1);
    
    atChargingStation=sum(u(i+1,:)==chargingStations(:,1));
    
    % update current statuses
    s(1,:)=logical(atChargingStation.*(d(i+1,:)==0));
    s(2,:)=logical(~atChargingStation.*(d(i+1,:)==0));
    s(3,:)=(d(i+1,:)>0).*s(3,:);
    s(4,:)=(d(i+1,:)>0).*s(4,:);
    s(5,:)=(d(i+1,:)>0).*s(5,:);
    
    % record current status
    status(i,:)=(1:5)*s;
    
    % record time
    S.lasttime=cputime;
    
    % update cputime of this step
	S.comptime(i+1)=cputime;
    
end


%% final results

% vehicle related
Sim.u=uint16(u); % final destination of vehicles (station) [tsim x m]
Sim.q=single(q); % state of charge 
Sim.e=datacompactor(e*P.Tech.battery*60/P.Sim.e);
Sim.ef=datacompactor(ef*P.Tech.battery*60/P.Sim.e);
Sim.status=status;

% passenger related
Sim.waiting=sparse(waiting); % waiting times
Sim.dropped=sparse(dropped); % dropped requests
Sim.chosenmode=chosenmode; % chosen mode
Sim.waitingestimated=sparse(waitingestimated); % estimated waiting time (only mode choice)
Sim.modalshare=sum(chosenmode)/r;

% general info
Sim.relodist=relodist*P.Sim.e; % relocation minutes
Sim.relodistkm=relodistkm; % relocation km
Sim.tripdist=tripdist*P.Sim.e; % trip minutes
Sim.tripdistkm=tripdistkm; % trip minutes
Sim.emissions=(sum(Sim.e/60*P.Sim.e,2)')*co2(1:tsim)/10^6; % emissions [ton]

% pricing info
Sim.revenues=sum((offeredprices-P.Pricing.movingcostkm.*tripDistancesKm(1:r)).*chosenmode.*(1-dropped));
Sim.relocationcosts=sum(relodistkm)*P.Pricing.movingcostkm;
Sim.offeredprices=datacompactor(offeredprices);
Sim.tariff=datacompactor(tariff);
Sim.surcharge=datacompactor(surcharge);

% Internals struct
Internals.b=b;
Internals.d=uint8(d);
Internals.zmacro=zmacro;


%% create Res struct and save results

% parameters of simulation
Params.cumulativeTripArrivals=cumulativeTripArrivals;
Params.Tr=uint8(Tr);
Params.elep=elep;
Params.co2=co2;
Params.tsim=tsim;

% Res struct generation
Res.Params=Params;
Res.Trips=Trips;
Res.Sim=Sim;
Res.Internals=Internals;
Res.CPUtimes=S;
Res.cputime=cputime-S.starttime;
Res.cost=(sum(Sim.e/60/1000*P.Sim.e,2)')*elep(1:tsim)+Sim.emissions*P.carbonprice;
Res.dropped=sum(dropped)/r;
Res.peakwait=max(waiting);
Res.avgwait=mean(waiting);



