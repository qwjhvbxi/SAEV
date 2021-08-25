%% [Res]=MAINSIM(P,dispiter)
% Main simulation.
% Vehicles start at beginning of time step, arrive at end of time step.
% 
% Problem: pricing cannot be run with charging because of how EMD is
% calculated right now
% Problem: OD-based pricing cannot run with aggregate prediction
% TODO: - T should be a matrix not struct
%       - can avoid using Par, instead just use P?
% 
% See also: main

function [Res]=mainsim(P,dispiter)

DataFolder=getdatafolder();

%% load external files: scenario, trips 

% load distance matrix
[T,D,clusters,clusterCenters,chargingStations]=getscenario(P.scenario);

% load trips
[A,Atimes,cumulativeTripArrivals,~]=gettrips(P);

% load electricity prices and carbon emissions
[elepMinute,co2Minute,~]=readexternalfile([DataFolder 'grid/' P.gridfile '.csv'],P.gridday,true);


%% setup clustering

n=size(clusters,1);             % number of nodes 
As=clusters(A(:,1:2));          % OD at clusters
nc=length(clusterCenters);      % number of clusters
ncs=length(chargingStations);   % number of CS


%% load predictions

% Aforecast is a t x n^2 matrix with probabilities of each OD pair at each
% time interval (needed only for pricing optimization/modechoice)
[fo,fd,dkod,Aforecast]=loadpredictions(P,As,Atimes);
% [fo,fd,dkod,dkemd,Aforecast]=loadpredictions(P,As,Atimes,Ts);


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
queue=zeros(100,1);         % temporary variable to store queued arrivals
g=zeros(1,P.m);             % current idle time
s=false(5,P.m);             % current status: 
%     1 charging
%     2 idle
%     3 moving
%     4 moving for relocation
%     5 moving to charging station

% results variables
waitingestimated=zeros(r,1);    % estimated minutes to wait for each request
waiting=zeros(r,1);             % minutes waited for each request
dropped=false(r,1);             % request is dropped?
chosenmode=false(r,1);          % which mode is chosen?
offeredprices=zeros(r,1);       % price offered to each passenger
modeutilities=zeros(r,2);       % utility of each mode from passenger [saev, alternative]
status=zeros(tsim,P.m);         % vehicle status
relodist=zeros(tsim,1);         % distances of relocation (at moment of decision)
relodistkm=zeros(tsim,1);       % distances of relocation in km (at moment of decision)
tripdist=zeros(tsim,1);         % distances of trips (at moment of acceptance)
tripdistkm=zeros(tsim,1);       % distances of trips in km (at moment of acceptance)
% pooling=zeros(r,1);             % pool ID of each user (if ride shared)


%% setup internal parameters

Par=struct('D',D,'Epsilon',P.Sim.e,'minsoc',P.Operations.minsoc,'maxsoc',P.Operations.maxsoc,'modechoice',P.modechoice,...
    'battery',P.Tech.battery,'maxwait',P.Operations.maxwait,'VOT',P.Pricing.VOT,'traveltimecost',P.Pricing.traveltimecost,...
    'LimitFCR',0,'chargepenalty',1,'v2gminsoc',P.Operations.v2gminsoc,'efficiency',P.Tech.efficiency,'SPlength',1,...
    'cyclingcost',P.Tech.cyclingcost,'carbonprice',P.carbonprice,'v2g',P.Operations.v2g,'fastchargesoc',0,'slowchargeratio',1,...
    'refillmaxsoc',0,'aggregateratio',1,'chargekw',P.Tech.chargekw,'consumption',P.Tech.consumption,'cssize',chargingStations(:,2));

% distances
if ~isstruct(T)
    Tr=max(1,round(T/P.Sim.e));% distance matrix in steps
    Tr(1:n+1:end)=0;          % no distance between same node
    Trs=Tr(clusterCenters,clusterCenters);
    Ts=T(clusterCenters,clusterCenters);
    Par.Tr=Tr;
    [~,closestCS]=min(Tr(:,chargingStations(:,1)),[],2); % closest charging station to each node
else
    Ts=T;
    for i=1:length(T)
        Ts(i).traveltime=T(i).traveltime(clusterCenters,clusterCenters);
        Ts(i).hour=T(i).hour;
    end
end
tripDistances=nan(r,1);
tripDistancesKm=D(sub2ind(size(D),A(:,1),A(:,2)));


%% setup relocation module

if isfield(P,'Relocation') && ~isempty(P.Relocation)
    ts=round(P.Relocation.ts/P.Sim.e);
    tr=round(P.Relocation.tr/P.Sim.e);
    tx=round(P.Relocation.tx/P.Sim.e);
    b=zeros(ceil(tsim/tx),nc,'double');  % imbalance
    autoRelocation=1;
else
    b=NaN;
    autoRelocation=0;
end


%% setup charging module

dynamicCharging=false;
Beta=0;
setPoints=[];

if isfield(P,'Charging') && ~isempty(P.Charging)
    
    if strcmp(P.Charging,'night')
        
        % night charging
        limitHour=5;        % charging hour limit
        Par.refillmaxsoc=0.6;   % min. soc to trigger charging during day
        zmacro=[ [1;0;1;0]*ones(1,60*limitHour/P.Sim.e) , [0;0;1;0]*ones(1,60*(24-limitHour)/P.Sim.e) ];
        
%         TODO: both night and non dynamic charge as much as possible, but for
%         night, these parameters change after a certain time:
%         Par.slowchargeratio=0;
%         Par.fastchargesoc=0.6;
        
    else
    
        dynamicCharging=true;
        
        % inputs
        Beta=P.Charging.beta;
        chargingHorizon=round(P.Charging.mthor/Beta);
        
        % variables
        etsim=floor(1440/Beta); % number of charging decisions
        melep=average2(elepMinute,Beta);
        mco2=average2(co2Minute,Beta);
        
        % energy layer variable: static values
        Par.Beta=Beta;
        Par.chargingHorizon=chargingHorizon;
        Par.minsocfleet=Par.minsoc+P.Charging.extrasoc;
        
        % matrix of optimal control variables for energy layer
        zmacro=zeros(4,etsim+chargingHorizon); 
        
    end

else
    
    % charge as much as possible
    zmacro=[1;0;1;0]*ones(1,tsim);
    
end


%% setup frequency control reserve module

if isfield(P,'Charging') && isfield(P,'FCR') && ~isempty(P.FCR) 
    
    [fraw,~,fresolution]=readexternalfile([DataFolder 'grid/' P.FCR.filename '.csv'],P.gridday,false);
    f=average2(fraw,P.Sim.e/fresolution);
    
    Par.aggregateratio=P.FCR.aggregatechargeratio; % charge rate at aggregate level optimization
    Par.LimitFCR=ceil(P.FCR.contracted*1000/P.Tech.chargekw);
    Par.fcrcontracted=P.FCR.contracted;
    Par.fcrlimits=P.FCR.limits;
    Par.fastchargesoc=P.FCR.fastchargesoc;
    Par.slowchargeratio=P.FCR.slowchargeratio;
    
    Par.SPlength=Beta;
    
end


%% generate aggregate trip statistics

emdFileName=[DataFolder 'temp/emd-' P.tripfolder '-' num2str(P.tripday) '-' num2str(Beta) '.mat'];
if exist(emdFileName,'file')
    load(emdFileName,'dkemd','dkod');
else
    dkod=computetraveltime(A,Atimes,T,Beta); % TODO: should be on prediction
    dkemd=computeemd(fo,fd,Ts,Beta); 
    save(emdFileName,'dkemd','dkod');
end


%% setup pricing module 

addpath functions/pricing

% add info to Pricing struct
P.Pricing.relocation=autoRelocation;
P.Pricing.c=D(clusterCenters,clusterCenters);

% initializations
perDistanceTariff=ones(nc,nc).*P.Pricing.basetariffkm;  % matrix of fares
surchargeMat=zeros(nc,nc);                              % matrix of surcharges per stations

if isfield(P.Pricing,'alternativecost') && ~isempty(P.Pricing.alternativecost)
    Aaltp=P.Pricing.alternativecost;
elseif isfield(P.Pricing,'alternativecostfile') && ~isempty(P.Pricing.alternativecostfile)
    load([DataFolder 'trips/' P.tripfolder '/' P.Pricing.alternativecostfile],'alternativecost');
    Aaltp=alternativecost{P.tripday};
else
    % alternative price for each user calculated with trip distances in km
    Aaltp=P.Pricing.alternativecostkm*tripDistancesKm; 
end

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

% initial state of charge
q(1,:)=P.Operations.initialsoc.*ones(1,P.m);      

% initial position of vehicles
if isfield(P.Operations,'uinit')
    u(1,:)=P.Operations.uinit;
else
    u(1,:)=chargingStations(randi(nc,1,P.m),1);
end

% initial delay
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
S.comptime=[S.starttime;zeros(tsim,1)];
S.clocktime=zeros(tsim,1);


%% start of iterations

for i=1:tsim
    
    
	%% display progress
    
    startclock=tic;
    displayprogress(i,tsim,dispiter,sum(S.clocktime),40)
    ui=double(u(i,:));
    di=double(d(i,:));
    relodist(i)=0;
    relodistkm(i)=0;
    
    
    %% calculate current distance matrix
    
    if isstruct(T)
        [thisT]=gettraveltimenow(T,i*Par.Epsilon);
        Ts=thisT(clusterCenters,clusterCenters);
        Tr=max(1,round(thisT/P.Sim.e));% distance matrix in steps
        Trs=Tr(clusterCenters,clusterCenters);
        Par.Tr=Tr;
        [~,closestCS]=min(Tr(:,chargingStations(:,1)),[],2); % closest charging station to each node
    end
    
    
    %% move idle vehicles back to charging stations
    
    if nc<n
        idleTime=g.*(1-atChargingStation);
        IdleReached=(idleTime>=P.Operations.maxidle/P.Sim.e);
        ui(IdleReached)=chargingStations(closestCS(ui(IdleReached)),1);
        relodistCS=Tr(sub2ind(size(Tr),u(i,IdleReached),ui(IdleReached)));
        relodistCSkm=D(sub2ind(size(D),u(i,IdleReached),ui(IdleReached)));
        di(IdleReached)=relodistCS;
        s(2,IdleReached)=0;
        s(5,IdleReached)=1;
        relodist(i)=relodist(i)+sum(relodistCS);
        relodistkm(i)=relodistkm(i)+sum(relodistCSkm);
    end
    
    
    %% mode choice and pricing optimization

    if P.modechoice
        
        % launch pricing optimization
        if P.Pricing.dynamic && mod(i-1,tp)==0
            
            % TODO / PRICING: reorganize for imperfect prediction!
            %       price of alternative option should be dependent on OD, not single
            %       passenger (the one seen by optimization). Specific cost only for
            %       mode choice module
            
            % TODO / PRICING: at the end should be only OD pricing, no matter how it's
            % found! (OD based or node based)
            % TODO / PRICING: mintariff in pricing optimization module
            
            % current pricing number
            kp=ceil(i/tp);
            
            P.Pricing.T=Ts;
            
            % expected trips
            tripsExpected=cumulativeTripArrivals(i)+1:cumulativeTripArrivals(min(length(cumulativeTripArrivals),i+tpH+1));
            AsExpected=As(tripsExpected,:);     % current ODs
            alternativeCosts=Aaltp(tripsExpected);  % current costs for alternative mode 
            [a,Ib,~]=unique(AsExpected,'rows','stable');
            alternativeCostsMat=sparse(a(:,1),a(:,2),alternativeCosts(Ib),nc,nc);

            % launch pricing optimization
            [perDistanceTariff,surchargeNodes]=pricingmodule(P.Pricing,AsExpected,alternativeCostsMat,clusters(ui));

            tariff(:,kp)=perDistanceTariff(:);
            surcharge(:,kp)=surchargeNodes;

            surchargeMat=surchargeNodes(1:nc)+surchargeNodes(nc+1:2*nc)';
            
            % odpricing=perDistanceTariff(:)+surchargeNodes(1:nc)+surchargeNodes(nc+1:2*nc)';
            
        end

        % adjust relocation predictions with mode choice 
        if autoRelocation && mod(i-1,tx)==0
            
            % expected trips
            tripsExpected=cumulativeTripArrivals(i)+1:cumulativeTripArrivals(min(length(cumulativeTripArrivals),i+ts+tr));
            tripsExpected=unique(min(tripsExpected,length(Aaltp))); % remove trips not covered by pricing
            AsExpected=As(tripsExpected,:);     % current ODs
            alternativeCosts=Aaltp(tripsExpected);  % current costs for alternative mode 
            [a,Ib,~]=unique(AsExpected,'rows','stable');
            alternativeCostsMat=sparse(a(:,1),a(:,2),alternativeCosts(Ib),nc,nc);
        
            option1=exp(-max(P.Pricing.mintariff,perDistanceTariff.*P.Pricing.c)-surchargeMat);
            multiplier=option1./(option1+exp(-alternativeCostsMat));

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
        
        % available vehicles
        available=sum(s(1:2,:))';
        
        % Vin: vehicles information in the form: [station delay soc connected relocating]
        Vin=[clusters(ui) , Par.Epsilon*di' , available.*q(i,:)' , s(1,:)' , logical(s(4,:)+s(5,:))' ];
        
        % Nin: passenger information in the form: [number waiting, arrivals expected, departures expected]
        Nin=[dw , round(sum(fd(i:i+ts,:)))' , round(sum(fo(i:i+ts+tr,:)))'];
        
        Par.Trs=Ts;
        Par.limite=P.Relocation.ts;
        
        [Vout,bkt]=relocationmodule(Vin,Nin,Par);
        
        % update vehicles position
        used=logical(Vout(:,2));
        relodestinations=clusterCenters(Vout(used,1));
        relodisti=Tr(sub2ind(size(Tr),ui(used),relodestinations'));
        relodistkmi=D(sub2ind(size(D),ui(used),relodestinations'));
        
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
        pp=max(P.Pricing.mintariff , perDistanceTariff(selectorClusters).*tripDistancesKm(trips)); % +surchargeMat(selectorClusters); % TODO: fix!
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
        modeutilities(trips,:)=Bout(:,5:6);
        
        queue=trips(queuei(queuei>0));
        
    end
    
    
    %% charging optimization
    
    if dynamicCharging
        
        if rem(i,Beta/P.Sim.e)==1

            % index of energy layer
            t=(i-1)/(Beta/P.Sim.e)+1;
            
            dktripnow=dkod(t:t+chargingHorizon-1)+dkemd(t:t+chargingHorizon-1);    % minutes spent traveling during this horizon
            melepnow=melep(t:t+chargingHorizon-1)/1000; % convert to [$/kWh]
            mco2now=mco2(t:t+chargingHorizon-1); % [g/kWh]

            currentsp=chargingmodule(Par,q(i,:),dktripnow,melepnow,mco2now);
            
            zmacro(:,t)=currentsp;

        end
        
    else
        
        t=i;
        
    end
    
    
    %% simulation variables update
    
    mat1=(ui==chargingStations(:,1));
    atChargingStation=sum(mat1);
    whichcs=(1:ncs)*mat1;
    
%     Par.slowchargeratio=P.Charging.slowchargeratio(min(ceil(i/60*P.Sim.e),length(P.Charging.slowchargeratio)));
    
    if 1 %advancedCharging
        
        % advanced charging algorithm
        
        if rem(i-1,Par.SPlength/P.Sim.e)==0
        
            setPoints=setpointfleet(Par,q(i,:),s(1,:),whichcs,zmacro(1:2,t));
            
        end
        
        [ei,efi]=setpointvehicle(Par,q(i,:),s(1,:),whichcs,setPoints,f(i));
        
    else
    
        % basic legacy charging algorithm: does not consider station sizes or frequency reserve
        
        ei=simplecharging(Par,q(i,:),s(1,:),zmacro(1:3,t));
        efi=0;
    
    end
    
    e(i,:)=ei;
    ef(i,:)=efi;
    
    % update SOC 
    q(i+1,:)=q(i,:) + max(0,e(i,:)+ef(i,:)) + min(0,e(i,:)+ef(i,:))/P.Tech.efficiency  -(di>0).*(P.Tech.consumption/P.Tech.battery*P.Sim.e);

    % update position
    u(i+1,:)=ui;
    
    % update idle time
    g=(di==0).*(g+1);
    
    % update delay
    d(i+1,:)=max(0,di-1);
    
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
    
    % update cputime and clocktime of this step
	S.comptime(i+1)=cputime;
    S.clocktime(i)=toc(startclock);
    
end


%% final results

% passengers still in queue should be dropped
dropped(queue(queue>0))=1;

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
Sim.modeutilities=modeutilities;

% general info
Sim.relodist=relodist*P.Sim.e; % relocation minutes
Sim.relodistkm=relodistkm; % relocation km
Sim.tripdist=tripdist*P.Sim.e; % trip minutes
Sim.tripdistkm=tripdistkm; % trip km
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

% stats
Stats.vehicletrips=sum(chosenmode.*(1-dropped))/P.m; % trips per vehicle per day
Stats.vkt=(sum(tripdistkm)+sum(relodistkm))/P.m; % km driven per vehicle per day
Stats.evktshare=sum(relodistkm)/(sum(tripdistkm)+sum(relodistkm)); % share of empty km driven
Stats.runtime=(sum(tripdist)+sum(relodist))/60/P.m; % average vehicle use per day (hours)
Stats.avgtriplength=sum(tripdistkm)/sum(chosenmode.*(1-dropped));
Stats.avgrevenuevehicle=sum(Sim.revenues)/P.m;
Stats.dropped=sum(dropped)/r;
Stats.peakwait=max(waiting);
Stats.avgwait=mean(waiting);
Stats.chargingcost=(sum(Sim.e/60/1000*P.Sim.e,2)')*elep(1:tsim)+Sim.emissions*P.carbonprice;
Stats.modalshare=sum(chosenmode)/r;
Stats.cputime=cputime-S.starttime;
Stats.clocktime=sum(S.clocktime);
% Stats.cycles= % average battery cycles per vehicle per day % TODO: how to consider different SOC at start and end?


%% create Res struct and save results

% parameters of simulation
Params.cumulativeTripArrivals=cumulativeTripArrivals;
Params.Tr=uint8(Tr);
Params.elep=elep;
Params.co2=co2;
Params.tsim=tsim;

% Res struct generation
Res.Params=Params;
Res.Trips=struct('dkod',dkod,'dkemd',dkemd,'dktrip',dkod+dkemd);
Res.Sim=Sim;
Res.Stats=Stats;
Res.Internals=Internals;
Res.CPUtimes=S;
Res.cputime=cputime-S.starttime;
Res.cost=(sum(Sim.e/60/1000*P.Sim.e,2)')*elep(1:tsim)+Sim.emissions*P.carbonprice;
Res.dropped=sum(dropped)/r;
Res.peakwait=max(waiting);
Res.avgwait=mean(waiting);



