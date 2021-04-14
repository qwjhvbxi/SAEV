% [Res]=GENERALC(P[,extsave,dispiter])
% Run SAEV simulation and relocation/charge oprimization.
% Vehicles start at beginning of time step, arrive at end of time step.
% 
% u: destination or position at beginning of time step
% q: SOC at beginning of time step
% e: charging energy exchanged during time step
% ef: FCR energy exchanged during time step
% b: imbalance from relocation module
% d: delay at beginning of time step
% sc: status during time step: connected/not connected
% sm: status during time step: moving/not moving
% 
% see also CPAR

function [Res]=generalC(P,extsave,dispiter)

%% initializations and input check

addpath functions utilities 
DataFolder=setDataFolder();
if nargin<3
    dispiter=1;
    if nargin<2
        extsave=1;
    end
end


%% generate unique Hash for savefile

Hash=DataHash(P);
simname=[DataFolder 'out_saev/simulations/' Hash '.mat'];
if extsave<2
    if extsave>=0 && exist(simname,'file')
        % if the file already exists, just return the previous result
        load(simname,'Res');
        return
    end
end

if isfield(P,'trlayeralg') && strcmp(P.trlayeralg,'opti')
    Res=generalOpti(P,extsave,dispiter);
    return
end


%% load external files: scenario, trips 

% load distance matrix
load([DataFolder 'scenarios/' P.scenario '.mat'],'T','Clusters','chargingStations');

% load trips
% TODO: can add secondary trip file (real vs expected/forecasted)
[A,Atimes,AbuckC,~]=loadTrips(P);
AbuckC=AbuckC(1:P.Sim.e:end);

if ~P.Sim.mpcpredict
    temp1=strfind(P.tripfolder,'_');
    Pb.tripfolder=P.tripfolder(1:temp1(end)-1);
    Pb.tripday=P.tripday;
    Pb.ratio=str2double(P.tripfolder(temp1(end)+1:end)); % TODO: change!!
    [As2,~,AbuckC2,~]=loadTrips(Pb);
    AbuckC2=AbuckC2(1:P.Sim.e:end);
end
    

% load electricity prices and carbon emissions
[elepMinute,co2Minute,~]=ReadExtFile([DataFolder 'grid/' P.gridfile '.csv'],P.gridday,true);


%% parameters of simulation

% parameters
n=size(T,1);              % number of nodes
r=AbuckC(1441);              % number of requests
tsim=1440/P.Sim.e;            % number of time steps
Tr=max(1,round(T/P.Sim.e));   % distance matrix in steps
Tr(1:n+1:end)=0;          % no distance between same node
aggregateratio=1;         % charge rate at aggregate level optimization
ac=P.Tech.chargekw/P.Tech.battery/60*P.Sim.e;    % charge rate per time step (normalized)
ad=P.Tech.consumption/P.Tech.battery*P.Sim.e;    % discharge rate per time step (normalized)
TripDistances=Tr(sub2ind(size(Tr),A(:,1),A(:,2)))*P.Sim.e; % trip distances in minutes

% electricity and emissions profiles
elep=average2(elepMinute,P.Sim.e);
co2=average2(co2Minute,P.Sim.e);

% check optional info
if isfield(P,'Pricing') && ~isempty(P.Pricing)
    Pricing=P.Pricing;
else
    Pricing=struct('relocationcost',0,'basetariff',0,'altp',0,'VOT',0,'pricingwaiting',1,'alternative',0,'dynamic',0);
end

% main simulation variables
q=zeros(tsim,P.m); % SOC
u=zeros(tsim,P.m,'uint16'); % vehicles in charging stations
d=zeros(tsim,P.m,'uint16'); % delay
e=zeros(tsim,P.m); % charging
ef=zeros(tsim,P.m); % FCR charging

% working variables
g=zeros(1,P.m); % current idle time
s=false(3,P.m); % current status: [relocating, connected, moving to charging station]
queue=zeros(100,1);         % temporary variable to store queued arrivals
DayCharge=0; % wether to charge during the day for unscheduled charging

% results variables
waitingestimated=zeros(r,1);    % estimated minutes to wait for each request
waiting=zeros(r,1);             % minutes waited for each request
dropped=false(r,1);             % request is dropped?
chosenmode=false(r,1);          % which mode is chosen?
% pooling=zeros(r,1);             % pool ID of each user (if ride shared)
offeredprices=zeros(r,1);       % price offered to each passenger
sc=false(tsim,P.m);             % connected to charging station status
sm=false(tsim,P.m);             % moving status
relodist=zeros(tsim,1);         % distances of relocation (at moment of decision)
tripdist=zeros(tsim,1);         % distances of trips (at moment of acceptance)


%% setup clustering

% override loaded info if passed from P struct
if isfield(P,'clusters')
    chargingStations=P.chargingStations;
    Clusters=P.clusters;
end

if exist('Clusters','var')
    As=Clusters(A);
    Trs=Tr(chargingStations,chargingStations);
    nc=length(chargingStations);             % number of clusters
else
    chargingStations=(1:n)';
    Clusters=(1:n)';
    As=A;
    Trs=Tr;
    nc=n;
end


%% setup relocation module

if isfield(P,'Relocation') && ~isempty(P.Relocation)
    ts=round(P.Relocation.ts/P.Sim.e);
    tr=round(P.Relocation.tr/P.Sim.e);
    tx=round(P.Relocation.tx/P.Sim.e);
    bmin=P.Relocation.bmin;
    b=zeros(ceil(tsim/tx),nc,'double');             % imbalance
    AutoRelo=1;
else
    AutoRelo=0;
end


%% setup charging module

DynamicCharging=false;
refillmaxsoc=0;
Trips=[];

if isfield(P,'Charging') && ~isempty(P.Charging)
    
    if strcmp(P.Charging,'night')
        
        LimitHour=5;
        
        % night charging
        zmacro=[ [1;0;1;0]*ones(1,60*LimitHour/P.Sim.e) , [0;0;1;0]*ones(1,60*(24-LimitHour)/P.Sim.e) ];
        
        % refill low soc
        refillmaxsoc=0.6;
        
    else
    
        DynamicCharging=true;
        
        % inputs
        Beta=P.Charging.beta;
        mthor=round(P.Charging.mthor/Beta);
        
        % variables
        etsim=floor(1440/Beta); % number of charging decisions
        melep=average2(elepMinute,Beta);
        mco2=average2(co2Minute,Beta);

        % generate aggregate trip statistics
        % EMDFileName=[P.tripfolder '-' num2str(P.tripday)];
        % Trips=generateEMD(A,Atimes,T,Beta,EMDFileName);
            
        EMDFileName=[P.tripfolder '-' num2str(P.tripday) '-50'];
        Trips=generateEMD(A,Atimes,T,Beta,EMDFileName,chargingStations,Clusters);
        
        % energy layer variable: static values
        E.v2g=P.Operations.v2g; % use V2G?
        E.eta=P.Tech.efficiency;       % roundtrip (discharge) efficiency
        E.selling=1;            % can sell to the grid?
        % TODO: socboost should be calculated by optimization as function of inputs
        E.socboost=1e4;         % soft constraint for final soc 
        E.T=mthor;% number of time steps in energy layer
        E.cyclingcost=P.Tech.cyclingcost;                       % battery cycling cost [$/kWh]
        E.storagemax=P.Tech.battery*P.m*P.Operations.maxsoc;    % max total energy in batteries [kWh]
        E.carbonprice=P.carbonprice;                            % carbon price [$ per kg]

        % matrix of optimal control variables for energy layer
        zmacro=zeros(4,etsim+mthor); 
        
    end

else 
    
    % charge as much as possible
    zmacro=[1;0;1;0]*ones(1,tsim);
    
end


%% setup frequency control reserve module

if isfield(P,'Charging') && isfield(P,'FCR') && ~isempty(P.FCR) 
    
    addpath functions/FCR
    
    FCR=true;
    if isfield(P.FCR,'aggregatechargeratio')
        aggregateratio=P.FCR.aggregatechargeratio;
    end
    
    % TODO: use CSV as inputs
    [fraw,~,fresolution]=ReadExtFile([DataFolder 'grid/' P.FCR.filename],P.gridday,false);
    f=average2(fraw,P.Sim.e/fresolution);
%     load([DataFolder 'grid/' P.FCR.filename],'f');
%     ReshapeFactor=size(f,1)/1440*P.Sim.e;
%     f=average2(f(:,P.gridday),ReshapeFactor);
    
    LimitFCR=ceil(P.FCR.contracted*1000/P.Tech.chargekw);
    
    ParC.ac=ac;
    ParC.af=P.FCR.contracted*1000/P.Tech.battery/60*P.Sim.e;    % FCR rate per time step (normalized)
    ParC.H=P.Charging.beta/P.Sim.e;
    ParC.limits=P.FCR.limits;
    ParC.fastchargesoc=P.FCR.fastchargesoc;
    ParC.slowchargeratio=P.FCR.slowchargeratio;
    ParC.minsoc=P.Operations.minsoc;
    ParC.maxsoc=P.Operations.maxsoc;
    ParC.v2gminsoc=P.Operations.v2gminsoc;
    ParC.battery=P.Tech.battery;
    ParC.efficiency=P.Tech.efficiency;
    
else
    FCR=false;
    LimitFCR=0;
end


%% setup pricing module 

% add info to Pricing struct
Pricing.c=Trs*P.Sim.e;
Pricing.relocation=AutoRelo;

% initialize matrix of fare per minute
PerDistanceTariff=ones(nc,nc).*Pricing.basetariff;

% initialize surcharges per stations
Surcharges=zeros(nc,nc);

% price of alternative option 
if numel(Pricing.alternative)>1
    Aaltp=Pricing.alternative; % alternative price for each user is given as input 
else
    Aaltp=Pricing.alternative*TripDistances; % alternative price for each user
    altpmat=Pricing.alternative.*Pricing.c; % alternative price for each OD
%     Pricing.altp=Pricing.alternative.*Pricing.c; % alternative price for each OD
end

% initialize multiplier for relocation
Multiplier=1;

if Pricing.dynamic
    addpath functions/pricing
    tp=round(Pricing.tp/P.Sim.e);       % pricing interval
    tpH=round(Pricing.horizon/P.Sim.e); % pricing horizon
    nodebased=isfield(Pricing,'nodebased') && Pricing.nodebased;
    
    % initialization for dynamic pricing
    if ~nodebased
        % initialize matrix of real prices offered
        DynamicPricing=ones(nc,nc,ceil(tsim/tp)+1)*Pricing.basetariff;
    else
        % initialize matrix of surcharges
        DynamicPricing=zeros(ceil(tsim/tp)+1,nc*2);
    end
    
    % dynamic calculation is ignored if modechoice is not selected
    if ~P.modechoice
        warning('modechoice option is disabled, ignoring price optimization.')
    end
else 
    tp=1;
    DynamicPricing=NaN;
end

% function to calculate probability of acceptance given a certain price for each OD pair
ProbAcc=@(p,s,altp) exp(-p.*Pricing.c-s)./(exp(-p.*Pricing.c-s)+exp(-altp));


%% setup trip assignment module

Par=struct('Tr',Tr,'ad',ad,'e',P.Sim.e,'minsoc',P.Operations.minsoc,'modechoice',P.modechoice,...
    'maxwait',P.Operations.maxwait,'VOT',Pricing.VOT,'WaitingCostToggle',Pricing.pricingwaiting,'LimitFCR',LimitFCR,'chargepenalty',1);

% if DynamicCharging
%     Par.chargepenalty=1;
% end


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

AtChargingStation=sum(u(1,:)==chargingStations);
    
% update statuses
s(2,:)=logical(AtChargingStation.*(d(1,:)==0));


%% variables for progress display and display initializations

S.starttime=cputime;
S.lasttime=S.starttime;
S.comptime=[cputime;zeros(tsim,1)];


%% start of iterations

for i=1:tsim
    
    
	%% display progress
    
    displayState(i,tsim,dispiter,S.comptime(i)-S.comptime(1),40)
    ui=double(u(i,:));
    di=double(d(i,:));
    relodist(i)=0;
    
    
    %% move idle vehicles back to charging stations
    
    if nc<n
        IdleReached=(g.*(1-AtChargingStation)>=P.Operations.maxidle/P.Sim.e);
        ui(IdleReached)=chargingStations(Clusters(ui(IdleReached)));
        relodistCS=Tr(sub2ind(size(Tr),u(i,IdleReached),ui(IdleReached)));
        di(IdleReached)=relodistCS;
        s(3,IdleReached)=true;
        relodist(i)=relodist(i)+sum(relodistCS);
    end
    
    
    %% charging optimization
    
    if DynamicCharging 
        
        if rem(i,Beta/P.Sim.e)==1
        
            % index of energy layer
            t=(i-1)/(Beta/P.Sim.e)+1;

            actualminsoc=min(P.Operations.minsoc+P.Charging.extrasoc,mean(q(i,:))*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of energy layer
            dktripnow=Trips.dktrip(t:t+mthor-1);    % minutes spent traveling during this horizon

            E.maxchargeminute=P.Tech.chargekw/60*aggregateratio;    % energy exchangeable per minute per vehicle [kWh]
            E.storagemin=P.Tech.battery*P.m*actualminsoc; % kWh
            E.einit=sum(q(i,:))*P.Tech.battery;     % total initial energy [kWh]
            E.etrip=dktripnow*P.Tech.consumption;   % energy used per step [kWh] 
            E.dkav=max(0,P.m*Beta-dktripnow);         % minutes of availability of cars
            E.electricityprice=melep(t:t+mthor-1)/1000; % convert to [$/kWh]
            E.emissionsGridProfile=mco2(t:t+mthor-1); % [g/kWh]

            maxc=E.dkav*E.maxchargeminute; % max exchangeable energy per time step [kWh]

            ELayerResults=aevopti11(E);

            if ~isempty(ELayerResults)

                % make sure that there is no discharging when V2G is not allowed
                if ~P.Operations.v2g
                    ELayerResults.discharging(:)=0;
                end

                % zmacro: [relative charging, relative discharging, max charging allowed, energy required by trips]
                zmacro(:,t:t+mthor-1)=[ ELayerResults.charging , ELayerResults.discharging , maxc , E.etrip ]';
                zmacro(isnan(zmacro))=0;

            else

                % in case there is no feasible solution, charge as much as possible
                zmacro(:,t:t+mthor-1)=[ones(size(maxc,1),1),zeros(size(maxc,1),1),maxc,E.etrip]';

            end


            % record time
            S.lasttime=cputime;

        end
    
    else
        
        t=i;
        
    end
    
    v2gallowed=q(i,:)>P.Operations.v2gminsoc;
    extracharge=(q(i,:)<refillmaxsoc);
    chargevector=max(-1,min(1,(ones(1,P.m)*(zmacro(1,t)/zmacro(3,t))-v2gallowed*(zmacro(2,t)/zmacro(3,t))+extracharge)))*ac;
    
    
    %% pricing optimization

    % current pricing number
    kp=ceil(i/tp);

    if P.modechoice
        
        % if Pricing.dynamic>0 && dynamicpricing<3 && (i==1 || mod(i-(ts+tr+1),tp)==0)
        if Pricing.dynamic
            
            if mod(i-1,tp)==0
        
                % expected trips
                Selection0=AbuckC(i)+1:AbuckC(min(length(AbuckC),i+tpH));
                AsNow=As(Selection0,:);
                a_tp=sparse(AsNow(:,1),AsNow(:,2),1,nc,nc);%+q_t;

                % number of vehicles at each station (including vehicles directed there)
                Pricing.v=histc(Clusters(ui),1:nc);
                
                a_tp(1:nc+1:end)=0;
                Pricing.a=a_tp;
                
                if numel(P.Pricing.alternative)>1
                    altpNow=Aaltp(Selection0);
                    [a,Ib,~]=unique(AsNow,'rows','stable');
                    Pricing.altp=sparse(a(:,1),a(:,2),altpNow(Ib),nc,nc);
                else
                    Pricing.altp=altpmat;
                end
                
                if ~nodebased
            
                    [ODpricesNow,~,~]=NLPricing5(Pricing); % OD-pair-based pricing

                    DynamicPricing(:,:,kp)=ODpricesNow;

                    PerDistanceTariff=DynamicPricing(:,:,kp);
                    
                else 
                    
                    [NodeSurchargeNow,~,~]=NLPricingNodes(Pricing); % node-based pricing

                    DynamicPricing(kp,:)=NodeSurchargeNow';

                    Surcharges=DynamicPricing(kp,1:nc)+DynamicPricing(kp,nc+1:2*nc)';
                
                end
            end  
        else
            
            % expected trips
            Selection0=AbuckC(i)+1:AbuckC(i+1);
            AsNow=As(Selection0,:);
            if numel(P.Pricing.alternative)>1
                altpNow=Aaltp(Selection0);
                [a,Ib,~]=unique(AsNow,'rows','stable');
                Pricing.altp=sparse(a(:,1),a(:,2),altpNow(Ib),nc,nc);
            else
                Pricing.altp=altpmat;
            end
            
        end
        
        Multiplier=(ProbAcc(PerDistanceTariff,Surcharges,Pricing.altp));
        
    end
    

    %% relocation

    % if it's time for a relocation decision
    if AutoRelo && mod(i-1,tx)==0
        
        % current relocation number
        kt=(i-1)/tx+1;

        % number of waiting passenger at station
        if sum(queue)>0
            dw=histcounts(As(queue(queue>0),1),1:nc+1)';
        else 
            dw=zeros(nc,1);
        end
        
        if P.Sim.mpcpredict
        
            % expected trips
            Selection1=AbuckC(i)+1:AbuckC(min(length(AbuckC),i+ts));
            a_ts=(Multiplier.*sparse(As(Selection1,1),As(Selection1,2),1,nc,nc));

            Selection2=AbuckC(i)+1:AbuckC(min(length(AbuckC),i+ts+tr));
            a_to=(Multiplier.*sparse(As(Selection2,1),As(Selection2,2),1,nc,nc));
        
        else
            
            % stochastic prediction
            
            % expected trips
            Selection1=AbuckC2(i)+1:AbuckC2(min(length(AbuckC2),i+ts));
            a_ts=(Multiplier.*sparse(As2(Selection1,1),As2(Selection1,2),1,nc,nc))/Pb.ratio;

            Selection2=AbuckC2(i)+1:AbuckC2(min(length(AbuckC2),i+ts+tr));
            a_to=(Multiplier.*sparse(As2(Selection2,1),As2(Selection2,2),1,nc,nc))/Pb.ratio;
            
        end

        % Vin: vehicles information in the form: [station delay soc connected relocating]
        Vin=[Clusters(ui) , di' , q(i,:)' , s(2,:)' , logical(s(1,:)+s(3,:))' ];
        ParRel.ad=ad;
        ParRel.dw=dw; % number of passengers waiting at each station
        ParRel.a_ts=round(sum(a_ts))'; % expected arrivals between now and now+ts
        ParRel.a_to=round(sum(a_to,2)); % expected requests between now and now+ts+tr
        ParRel.Trs=Trs;
        ParRel.limite=P.Relocation.ts;
        ParRel.bmin=bmin;
        ParRel.LimitFCR=LimitFCR;
        ParRel.chargepenalty=Par.chargepenalty;
        
        [Vout,bkt]=Relocation(Vin,ParRel);
        
        % update vehicles position
        used=logical(Vout(:,2));
        relodestinations=chargingStations(Vout(used,1));
        relodisti=Tr(sub2ind(size(Tr),ui(used),relodestinations'));
        
        ui(used)=relodestinations; 
        di(used)=di(used)+relodisti;
        
        % update vehicles status (relocating vehicles cannot be relocated)
        s(1,used)=1;
        s(2:3,used)=0;
        
        % update results
        b(kt,:)=bkt;
        relodist(i)=relodist(i)+sum(relodisti);
        
    end


    %% trip assignment

    % generate trip requests for this time step
    trips=(AbuckC(i)+1:AbuckC(i+1))';

    % add previously queued requests and reset queue
    if sum(queue)>0
        trips=[queue(queue>0);trips];
        queue(:)=0;
    end
    
    if ~isempty(trips) 
    
        % TODO: fix pooling option 
        
        % calculate pricing    
        SelectorClusters=sub2ind(size(Trs),As(trips,1),As(trips,2));
        pp=PerDistanceTariff(SelectorClusters).*TripDistances(trips)+...
            Surcharges(SelectorClusters);
        alte=exp(-Aaltp(trips));
        
        % offered prices
        offeredprices(trips)=pp;
        
        % Vin: vehicles information in the form: [station delay soc connected]
        Vin=[ui' , di' , q(i,:)' , s(2,:)' ];
        
        % Bin: passengers info in the form: [O D waiting offeredprice utilityalternative]
        Bin=[A(trips,:) , waiting(trips) , pp , alte ];

        % tripAssignment (no clustering) or tripAssignment2 (clustering)
        if AutoRelo
            [Vout,Bout,tripdisti,relodistiPU,queuei]=tripAssignment2(Vin,Bin,Par);
        else
            [Vout,Bout,tripdisti,relodistiPU,queuei]=tripAssignment(Vin,Bin,Par);
        end
        
        % update vehicles positions
        ui=Vout(:,1)';
        di=Vout(:,2)';
        
        % update vehicles status
        used=logical(Vout(:,3));
        s(1:3,used)=0;
        
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
    
    if FCR  % setpoint based
        
        if rem(i,Beta/P.Sim.e)==1
            
            [SetPoints]=fleetSetpoint(ParC,q(i,:),s(2,:),zmacro(1:2,t));
            
        end
        
        [ei,efi]=vehicleSetpoint(ParC,q(i,:),s(2,:),f(i),SetPoints);
        
        e(i,:)=ei;
        ef(i,:)=efi;
           
    else
        
        % capacity based
        CapUp=s(2,:).*min(ac,P.Operations.maxsoc-q(i,:)); % charge
        CapDown=s(2,:).*min(ac,(q(i,:)-P.Operations.minsoc)*P.Tech.efficiency); % discharge
        
        e(i,:)=min(CapUp,max(0,chargevector))+max(-CapDown,min(0,chargevector));
        
    end

    % update SOC 
    q(i+1,:)=q(i,:)+max(0,e(i,:)+ef(i,:))+min(0,e(i,:)+ef(i,:))/P.Tech.efficiency-(di>0).*ad;

    % update position
    u(i+1,:)=ui;
    
    % update idle time
    g=(di==0).*(g+1);
    
    % update delay
    d(i+1,:)=max(0,di-1);
    
    AtChargingStation=sum(u(i+1,:)==chargingStations);
    
    % update historic statuses
    sm(i,:)=logical(di>0);
    sc(i,:)=logical(s(2,:));
    
    % update current statuses
    s(1,:)=(d(i+1,:)>0).*s(1,:);
    s(2,:)=logical(AtChargingStation.*(d(i+1,:)==0));
    s(3,:)=(d(i+1,:)>0).*s(3,:);
    
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
if FCR
    Sim.ef=single(ef*P.Tech.battery*60/P.Sim.e);
else
    Sim.ef=sparse(ef);
end

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
if isfield(P,'Pricing')
    Sim.revenues=sum((offeredprices-Pricing.relocationcost.*TripDistances).*chosenmode.*(1-dropped));
    Sim.relocationcosts=sum(relodist)*P.Sim.e*Pricing.relocationcost;
    Sim.offeredprices=offeredprices;
    Sim.prices=DynamicPricing;
end

% Internals struct
Internals.b=b;
Internals.sc=logical(sc);
Internals.sm=logical(sm);
Internals.d=uint8(d);
Internals.zmacro=zmacro;


%% create Res struct and save results

% total cpu time
elapsed=cputime-S.starttime;

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
Res.cputime=elapsed;
Res.cost=(sum(Sim.e/60/1000*P.Sim.e,2)')*elep(1:tsim)+Sim.emissions*P.carbonprice;
Res.dropped=sum(dropped)/r;
Res.peakwait=max(waiting);
Res.avgwait=mean(waiting);

% save results
if extsave>0
    save(simname,'Res','P');
end


%% end display

if dispiter<0
    fprintf('sim #%d ',-dispiter)
end
fprintf('successfully completed - avg soc: %0.2f - total time: %d:%0.2d - \n',mean(q(i,:)),floor(elapsed/60),round(rem(elapsed,60)));

end




