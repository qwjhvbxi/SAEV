%% [Res]=MAINSIM(P,dispiter)
% Main simulation.
% Vehicles start at beginning of time step, arrive at end of time step.
% 
% TODO: add charging station size
% 
% See also: main

function [Res]=mainsim(P,dispiter)

DataFolder=getdatafolder();

%% load external files: scenario, trips 

% load distance matrix
load([DataFolder 'scenarios/' P.scenario '.mat'],'T','Clusters','chargingStations');

% load trips
[A,Atimes,AbuckC,~]=gettrips(P);
AbuckC=AbuckC(1:P.Sim.e:end);

% load electricity prices and carbon emissions
[elepMinute,co2Minute,~]=readexternalfile([DataFolder 'grid/' P.gridfile '.csv'],P.gridday,true);


%% setup clustering

n=size(T,1);              % number of nodes
if exist('Clusters','var')
    clusters=Clusters;
else
    chargingStations=(1:n)';
    clusters=(1:n)';
end
As=clusters(A);
nc=length(chargingStations);             % number of clusters


%% initialize parameters of simulation

% parameters
r=AbuckC(1441);           % number of requests
tsim=1440/P.Sim.e;        % number of time steps
Tr=max(1,round(T/P.Sim.e));% distance matrix in steps
Tr(1:n+1:end)=0;          % no distance between same node
Ts=T(chargingStations,chargingStations);
Trs=Tr(chargingStations,chargingStations);
ac=P.Tech.chargekw/P.Tech.battery/60*P.Sim.e;    % charge rate per time step (normalized)
ad=P.Tech.consumption/P.Tech.battery*P.Sim.e;    % discharge rate per time step (normalized)
aggregateratio=1;         % charge rate at aggregate level optimization
tripDistances=Tr(sub2ind(size(Tr),A(:,1),A(:,2)))*P.Sim.e; % trip distances in minutes

% electricity and emissions profiles
elep=average2(elepMinute,P.Sim.e);
co2=average2(co2Minute,P.Sim.e);
f=sparse(tsim,1);

% check optional info
if isfield(P,'Pricing') && ~isempty(P.Pricing)
    Pricing=P.Pricing;
else
    Pricing=struct('relocationcost',0,'basetariff',0,'VOT',0,'pricingwaiting',1,'alternative',0,'dynamic',0);
end

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
tripdist=zeros(tsim,1);         % distances of trips (at moment of acceptance)
% pooling=zeros(r,1);             % pool ID of each user (if ride shared)


%% load predictions

% TODO: cleanup call to secondary trip file (real vs expected/forecasted)
if ~P.Sim.mpcpredict
    temp1=strfind(P.tripfolder,'_');
    Pb.tripfolder=P.tripfolder(1:temp1(end)-1);
    Pb.tripday=P.tripday;
    Pb.ratio=str2double(P.tripfolder(temp1(end)+1:end)); % TODO: change!!
    [As2,~,AbuckC2,~]=gettrips(Pb);
    AbuckC2=AbuckC2(1:P.Sim.e:end);
else
    Pb.ratio=1;
    As2=As;
    AbuckC2=AbuckC;
end 


%% setup trip assignment module

Par=struct('Tr',Tr,'ad',ad,'ac',ac,'e',P.Sim.e,'minsoc',P.Operations.minsoc,'maxsoc',P.Operations.maxsoc,'modechoice',P.modechoice,...
    'maxwait',P.Operations.maxwait,'VOT',Pricing.VOT,'WaitingCostToggle',Pricing.pricingwaiting,...
    'LimitFCR',0,'chargepenalty',1,'v2gminsoc',P.Operations.v2gminsoc,'efficiency',P.Tech.efficiency);

%% setup relocation module

if isfield(P,'Relocation') && ~isempty(P.Relocation)
    ts=round(P.Relocation.ts/P.Sim.e);
    tr=round(P.Relocation.tr/P.Sim.e);
    tx=round(P.Relocation.tx/P.Sim.e);
    bmin=P.Relocation.bmin;
    b=zeros(ceil(tsim/tx),nc,'double');             % imbalance
    autoRelocation=1;
else
    autoRelocation=0;
end


%% setup frequency control reserve module

Par.fcr=false;
Beta=0;
if isfield(P,'Charging') && isfield(P,'FCR') && ~isempty(P.FCR) 
    
    Par.fcr=true;
    if isfield(P.FCR,'aggregatechargeratio')
        aggregateratio=P.FCR.aggregatechargeratio;
    end
    
    % TODO: use CSV as inputs
    [fraw,~,fresolution]=readexternalfile([DataFolder 'grid/' P.FCR.filename],P.gridday,false);
    f=average2(fraw,P.Sim.e/fresolution);
    
    fcrLimit=ceil(P.FCR.contracted*1000/P.Tech.chargekw);
    Par.LimitFCR=fcrLimit;
    Par.af=P.FCR.contracted*1000/P.Tech.battery/60*P.Sim.e;    % FCR rate per time step (normalized)
    Par.H=P.Charging.beta/P.Sim.e;
    Par.limits=P.FCR.limits;
    Par.fastchargesoc=P.FCR.fastchargesoc;
    Par.slowchargeratio=P.FCR.slowchargeratio;
    Par.battery=P.Tech.battery;
    
end


%% setup charging module

dynamicCharging=false;
Par.refillmaxsoc=0;
Trips=[];

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
            dkod=computetraveltime(A,Atimes,T,Beta);
            dkemd=computeemd(As,Atimes,Ts,Beta); % TODO: should be on prediction
            dktrip=dkod+dkemd;
            save(emdFileName,'dkemd','dkod','dktrip');
        end
        Trips=struct('dkod',dkod,'dkemd',dkemd,'dktrip',dktrip);
        
        % energy layer variable: static values
        E.v2g=P.Operations.v2g; % use V2G?
        E.efficiency=P.Tech.efficiency;% roundtrip (discharge) efficiency
        E.socboost=1e4;         % soft constraint for final soc (TODO: change depending on inputs)
        E.T=chargingHorizon;              % number of time steps in energy layer
        E.cyclingcost=P.Tech.cyclingcost;                       % battery cycling cost [$/kWh]
        E.storagemax=P.Tech.battery*P.m*P.Operations.maxsoc;    % max total energy in batteries [kWh]
        E.carbonprice=P.carbonprice;                            % carbon price [$ per kg]
        E.maxchargeminute=P.Tech.chargekw/60*aggregateratio;    % energy exchangeable per minute per vehicle [kWh]
        
        % matrix of optimal control variables for energy layer
        zmacro=zeros(4,etsim+chargingHorizon); 
        
    end

else
    
    % charge as much as possible
    zmacro=[1;0;1;0]*ones(1,tsim);
    
end


%% setup pricing module 

% add info to Pricing struct
Pricing.n=nc;
Pricing.c=Trs*P.Sim.e;
Pricing.relocation=autoRelocation;

% price of alternative option 
if numel(Pricing.alternative)>1
    Aaltp=Pricing.alternative; % alternative price for each user is given as input 
else
    Aaltp=Pricing.alternative*tripDistances; % alternative price for each user
end

% initialize multiplier for relocation
multiplier=1;

if Pricing.dynamic
    
    addpath functions/pricing
    tp=round(Pricing.tp/P.Sim.e);       % pricing interval
    tpH=round(Pricing.horizon/P.Sim.e); % pricing horizon
    
    % initialization for dynamic pricing
    tariff=ones(nc^2,ceil(tsim/tp)+1)*Pricing.basetariff;
    surcharge=zeros(nc*2,ceil(tsim/tp)+1);
    
else
    
    tpH=0;
    tp=1;
    
end

% function to calculate probability of acceptance given a certain price for each OD pair
probAcc=@(p,s,altp) exp(-p.*Pricing.c-s)./(exp(-p.*Pricing.c-s)+exp(-altp));


%% initial states

q(1,:)=P.Operations.initialsoc.*ones(1,P.m);      % initial state of charge
if isfield(P.Operations,'uinit')
    u(1,:)=P.Operations.uinit;
else
    u(1,:)=chargingStations(randi(nc,1,P.m));                 % initial position of vehicles
end
if isfield(P.Operations,'dinit')
    d(1,:)=P.Operations.dinit;
end

% initial status
atChargingStation=sum(u(1,:)==chargingStations);
s(1,:)=logical(atChargingStation.*(d(1,:)==0));
s(2,:)=logical(~atChargingStation.*(d(1,:)==0));


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
    
    
    %% move idle vehicles back to charging stations
    
    if nc<n
        IdleReached=(g.*(1-atChargingStation)>=P.Operations.maxidle/P.Sim.e);
        ui(IdleReached)=chargingStations(clusters(ui(IdleReached)));
        relodistCS=Tr(sub2ind(size(Tr),u(i,IdleReached),ui(IdleReached)));
        di(IdleReached)=relodistCS;
        s(2,IdleReached)=0;
        s(5,IdleReached)=1;
        relodist(i)=relodist(i)+sum(relodistCS);
    end
    
    
    %% charging optimization
    
    if dynamicCharging 
        
        if rem(i,Beta/P.Sim.e)==1

            % index of energy layer
            t=(i-1)/(Beta/P.Sim.e)+1;

            actualminsoc=min(P.Operations.minsoc+P.Charging.extrasoc,mean(q(i,:))*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of energy layer
            dktripnow=dktrip(t:t+chargingHorizon-1);    % minutes spent traveling during this horizon

            E.storagemin=P.Tech.battery*P.m*actualminsoc; % kWh
            E.einit=sum(q(i,:))*P.Tech.battery;     % total initial energy [kWh]
            E.etrip=dktripnow*P.Tech.consumption;   % energy used per step [kWh] 
            E.dkav=max(0,P.m*Beta-dktripnow);         % minutes of availability of cars
            E.electricityprice=melep(t:t+chargingHorizon-1)/1000; % convert to [$/kWh]
            E.emissionsGridProfile=mco2(t:t+chargingHorizon-1); % [g/kWh]

            maxc=E.dkav*E.maxchargeminute; % max exchangeable energy per time step [kWh]

            ELayerResults=aevopti11(E);

            if ~isempty(ELayerResults)

                % zmacro: [relative charging, relative discharging, max charging allowed, energy required by trips]
                zmacro(:,t:t+chargingHorizon-1)=[ ELayerResults.charging , ELayerResults.discharging*logical(P.Operations.v2g) , maxc , E.etrip ]';
                zmacro(isnan(zmacro))=0;

            else

                % in case there is no feasible solution, charge as much as possible
                zmacro(:,t:t+chargingHorizon-1)=[ones(size(maxc,1),1),zeros(size(maxc,1),1),maxc,E.etrip]';

            end

            % record time
            S.lasttime=cputime;

        end
        
    else
        
        t=i;
        
    end
    
    
    %% pricing optimization

    if P.modechoice
        
        if ~Pricing.dynamic || mod(i-1,tp)==0
            
            % current pricing number
            kp=ceil(i/tp);
        
            % expected trips
            selection0=AbuckC(i)+1:AbuckC(min(length(AbuckC),i+tpH));
            AsNow=As(selection0,:);
            altpNow=Aaltp(selection0);
            
            [perDistanceTariff,surchargeNodes,altp]=pricingmodule(Pricing,AsNow,altpNow,clusters(ui));
            
            tariff(:,kp)=perDistanceTariff(:);
            surcharge(:,kp)=surchargeNodes;
            
            surchargeMat=surchargeNodes(1:nc)+surchargeNodes(nc+1:2*nc)';
            
        end
        
        multiplier=(probAcc(perDistanceTariff,surchargeMat,altp));
        
    end
    

    %% relocation

    % if it's time for a relocation decision
    if autoRelocation && mod(i-1,tx)==0
        
        % current relocation number
        kt=(i-1)/tx+1;

        % number of waiting passenger at station
        dw=histcounts(As(queue(queue>0),1),1:nc+1)';
            
        % expected trips
        selection1=AbuckC2(i)+1:AbuckC2(min(length(AbuckC2),i+ts));
        a_ts=(multiplier.*sparse(As2(selection1,1),As2(selection1,2),1,nc,nc))/Pb.ratio;

        Selection2=AbuckC2(i)+1:AbuckC2(min(length(AbuckC2),i+ts+tr));
        a_to=(multiplier.*sparse(As2(Selection2,1),As2(Selection2,2),1,nc,nc))/Pb.ratio;
            
        available=sum(s(1:2,:))';
        
        % Vin: vehicles information in the form: [station delay soc connected relocating]
        Vin=[clusters(ui) , di' , available.*q(i,:)' , s(1,:)' , logical(s(4,:)+s(5,:))' ];
        ParRel.dw=dw; % number of passengers waiting at each station
        ParRel.a_ts=round(sum(a_ts))'; % expected arrivals between now and now+ts
        ParRel.a_to=round(sum(a_to,2)); % expected requests between now and now+ts+tr
        ParRel.Trs=Trs;
        ParRel.limite=P.Relocation.ts;
        ParRel.bmin=bmin;
        ParRel.ad=ad;
        ParRel.LimitFCR=Par.LimitFCR;
        ParRel.chargepenalty=Par.chargepenalty;
        
        [Vout,bkt]=relocationmodule(Vin,ParRel);
        
        % update vehicles position
        used=logical(Vout(:,2));
        relodestinations=chargingStations(Vout(used,1));
        relodisti=Tr(sub2ind(size(Tr),ui(used),relodestinations'));
        
        ui(used)=relodestinations; 
        di(used)=di(used)+relodisti;
        
        % update vehicles status (relocating vehicles cannot be relocated)
        s(1:2,used)=0;
        s(4,used)=1;
        s(5,used)=0;
        
        % update results
        b(kt,:)=bkt;
        relodist(i)=relodist(i)+sum(relodisti);
        
    end


    %% trip assignment

    % generate trip requests for this time step
    trips=(AbuckC(i)+1:AbuckC(i+1))';

    % add previously queued requests and reset queue
    trips=[queue(queue>0);trips];
    queue(:)=0;
    
    if ~isempty(trips)
    
        % TODO: fix pooling option 
        
        % calculate pricing    
        selectorClusters=sub2ind(size(Trs),As(trips,1),As(trips,2));
        pp=perDistanceTariff(selectorClusters).*tripDistances(trips)+...
            surchargeMat(selectorClusters);
        alte=exp(-Aaltp(trips));
        
        % offered prices
        offeredprices(trips)=pp;
        
        % available vehicles depending on status
        available=sum(s(1:3,:))';
        
        % Vin: vehicles information in the form: [station delay soc connected]
        Vin=[ui' , di' , available.*q(i,:)' , s(1,:)' ];
        
        % Bin: passengers info in the form: [O D waiting offeredprice utilityalternative]
        Bin=[A(trips,:) , waiting(trips) , pp , alte ];

        if autoRelocation
            [Vout,Bout,tripdisti,relodistiPU,queuei]=tripassignmentsaev(Vin,Bin,Par);
        else
            [Vout,Bout,tripdisti,relodistiPU,queuei]=tripassignmentcarsharing(Vin,Bin,Par);
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
        relodist(i)=relodist(i)+relodistiPU;
        chosenmode(trips)=Bout(:,1);
        waiting(trips)=Bout(:,2);
        dropped(trips)=Bout(:,3);
        waitingestimated(trips)=waitingestimated(trips)+Bout(:,4);
        
        queue=trips(queuei(queuei>0));
        
    end
    
    
    %% simulation variables update
    
    [ei,efi]=chargingsetpoints(Par,q(i,:),s(1,:),zmacro(1:3,t),f(i),(rem(i,Beta/P.Sim.e)==1));
    e(i,:)=ei;
    ef(i,:)=efi;
    
    % update SOC 
    q(i+1,:)=q(i,:)+max(0,e(i,:)+ef(i,:))+min(0,e(i,:)+ef(i,:))/P.Tech.efficiency-(di>0).*ad;

    % update position
    u(i+1,:)=ui;
    
    % update idle time
    g=(di==0).*(g+1);
    
    % update delay
    d(i+1,:)=max(0,di-1);
    
    atChargingStation=sum(u(i+1,:)==chargingStations);
    
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
Sim.e=sparse(e/ac*P.Tech.chargekw);
Sim.ef=datacompactor(ef*P.Tech.battery*60/P.Sim.e);
Sim.status=status;

% passenger related
Sim.waiting=sparse(waiting); % waiting times
Sim.dropped=sparse(dropped); % dropped requests
Sim.chosenmode=chosenmode; % chosen mode
Sim.waitingestimated=sparse(waitingestimated); % estimated waiting time (only mode choice)

% general info
Sim.relodist=relodist*P.Sim.e; % relocation minutes
Sim.tripdist=tripdist*P.Sim.e; % trip minutes
Sim.emissions=(sum(Sim.e/60*P.Sim.e,2)')*co2(1:tsim)/10^6; % emissions [ton]

% pricing info
Sim.revenues=sum((offeredprices-Pricing.relocationcost.*tripDistances(1:r)).*chosenmode.*(1-dropped));
Sim.relocationcosts=sum(relodist)*P.Sim.e*Pricing.relocationcost;
Sim.offeredprices=datacompactor(offeredprices);
Sim.tariff=datacompactor(tariff);
Sim.surcharge=datacompactor(surcharge);

% Internals struct
Internals.b=b;
Internals.d=uint8(d);
Internals.zmacro=zmacro;


%% create Res struct and save results

% parameters of simulation
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



