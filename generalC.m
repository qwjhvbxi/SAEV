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

%% initializations

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

if strcmp(P.trlayeralg,'opti')
    Res=generalOpti(P,extsave,dispiter);
    return
end


%% parallel computing and input check

parcomp=is_in_parallel();

if isfield(P,'tripfolder')
    TripName=P.tripfolder;
else
    TripName=P.tripfile;
end


%% load external files: scenario, trips and energy

% load distance matrix
load([DataFolder 'scenarios/' P.scenario '.mat'],'T','Clusters','chargingStations');

% load trips
% TODO: can add secondary trip file (real vs expected/forecasted)
[A,Atimes,AbuckC,Distances]=loadTrips(P);
AbuckC=AbuckC(1:P.e:end);

% load electricity prices and carbon emissions
% x and y are electricity prices and carbon emissions, respectively. They
% are in matrix of shape N x days, with N the number of data points in a day
load([DataFolder 'eleprices/' P.gridfile '.mat'],'x','y');


%% parameters of simulation

% parameters
n=size(T,1);              % number of nodes
tsim=1440/P.e;            % number of time steps
etsim=floor(1440/P.beta); % number of charging decisions
Tr=max(1,round(T/P.e));   % distance matrix in steps
Tr(1:n+1:end)=0;
ac=P.Tech.chargekw/P.Tech.battery/60*P.e;    % charge rate per time step (normalized)
ad=P.Tech.consumption/P.Tech.battery*P.e;    % discharge rate per time step (normalized)
aggregateratio=1;         % charge rate at aggregate level optimization
mthor=round(P.EnergyLayer.mthor/P.beta);
if isfield(P.Tech,'efficiency')
    Efficiency=P.Tech.efficiency;
else
    Efficiency=1;
end
DayCharge=0;

% main variables
q=zeros(tsim,P.m,'double'); % SOC
u=zeros(tsim,P.m,'uint16'); % vehicles in charging stations
d=zeros(tsim,P.m,'uint16'); % delay
e=zeros(tsim,P.m,'double'); % charging
ef=zeros(tsim,P.m,'double'); % charging

g=zeros(1,P.m); % current idle time
s=false(3,P.m); % current status: [relocating, connected, moving to charging station]
sc=false(tsim,P.m); % connected to charging station status
sm=false(tsim,P.m); % moving status
queue=zeros(100,1);         % temporary variable to store queued arrivals

% results variables
waiting=zeros(length(A),1);  % minutes waited for each request
dropped=false(length(A),1);  % request is dropped?
chosenmode=false(length(A),1);% which mode is chosen?
pooling=zeros(length(A),1);  % pool ID of each user (if ride shared)
waitingestimated=zeros(length(A),1);  % estimated minutes to wait for each request
offeredprices=zeros(length(A),1);  % price offered to each passenger
relodist=zeros(tsim,1); % distances of relocation (at moment of decision)
tripdist=zeros(tsim,1); % distances of trips (at moment of acceptance)

% electicity price and emissions profiles
d1=P.gridday;
d2=rem(P.gridday,size(x,2))+1;
ReshapeFactor=tsim/size(x,1);
elep=repelem( [ x(:,d1);x(:,d2) ] , ReshapeFactor ,1); % electricity price profiles in $/MWh 
if exist('y','var') % carbon emissions profiles [g/kWh]
    co2=repelem( [ y(:,d1);y(:,d2) ] , ReshapeFactor ,1);
else
    co2=zeros(tsim*2,1);
end
melep=average2(elep,P.beta/P.e);
mco2=average2(co2,P.beta/P.e);
clear d1 d2 ReshapeFactor x y;

% trip distances
TripDistances=Tr(sub2ind(size(Tr),A(:,1),A(:,2)))*P.e; % minutes


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

if ~strcmp(P.trlayeralg,'no') && ~isempty(P.TransportLayer)
    ts=round(P.TransportLayer.ts/P.e);
    tr=round(P.TransportLayer.tr/P.e);
    tx=round(P.TransportLayer.tx/P.e);
    bmin=P.TransportLayer.bmin;
    b=zeros(ceil(tsim/tx),nc,'double');             % imbalance
    AutoRelo=1;
else
    ts=0;
    tr=0;
    b=NaN;
    AutoRelo=0;
end


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


%% setup charging module

if strcmp(P.enlayeralg,'aggregate') 
    
    % generate aggregate trip statistics
    EMDFileName=[TripName '-' num2str(P.tripday)];
    [Trips,~,~]=generateEMD(A,Atimes,T,etsim,EMDFileName);

    % append values for next day
    if isfield(P,'tripfolder')
        P2=P;
        P2.tripday=P.tripday+1;
        try 
            [A2,Atimes2,~,~]=loadTrips(P2);
        catch 
            P2.tripday=1;
            [A2,Atimes2,~,~]=loadTrips(P2);
        end
        EMDFileName=[TripName '-' num2str(P2.tripday)];
        [Trips2,~,~]=generateEMD(A2,Atimes2,T,etsim,EMDFileName);
        Trips.dktrip=[Trips.dktrip(1:48,:) ; Trips2.dktrip(1:48,:)];
    end

    % energy layer variable: static values
    E.v2g=P.Operations.v2g; % use V2G?
    E.eta=Efficiency;       % roundtrip (discharge) efficiency
    E.selling=1;            % can sell to the grid?
    % E.minfinalsoc=1;      % final SOC. This only works for optimization horizon of ~24h
    E.socboost=1e4;         % soft constraint for final soc
    E.T=mthor;% number of time steps in energy layer
    E.cyclingcost=P.Tech.cyclingcost;                       % battery cycling cost [$/kWh]
    E.storagemax=P.Tech.battery*P.m*P.Operations.maxsoc;    % max total energy in batteries [kWh]
    E.carbonprice=P.carbonprice;                            % carbon price [$ per kg]

    % matrix of optimal control variables for energy layer
    zmacro=zeros(4,etsim+mthor); 
    
else 
    
    Trips=[];
    zmacro=zeros(4,etsim);
    
end


%% frequency control reserve

if isfield(P,'FCR') && ~isempty(P.FCR)
    FCR=true;
    if isfield(P.FCR,'aggregatechargeratio')
        aggregateratio=P.FCR.aggregatechargeratio;
    end
    load([DataFolder 'grid/' P.FCR.filename],'f');
    ReshapeFactor=size(f,1)/1440*P.e;
    f=average2(f(:,P.gridday),ReshapeFactor);
    
    if isempty(P.FCR.contracted)
        
%         % dynamic variables
%         
%         E.maxchargeminute=P.Tech.chargekw/60*aggregateratio;    % energy exchangeable per minute per vehicle [kWh]
%         actualminsoc=min(P.Operations.minsoc+P.EnergyLayer.extrasoc,mean(q(1,:))*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of energy layer
%         E.storagemin=P.Tech.battery*P.m*actualminsoc; % kWh
%         dktripnow=Trips.dktrip(1:mthor);    % minutes spent traveling during this horizon
%         E.einit=sum(q(1,:))*P.Tech.battery;     % total initial energy [kWh]
%         E.etrip=dktripnow*P.Tech.consumption;   % energy used per step [kWh] 
%         E.dkav=max(0,P.m*P.beta-dktripnow);         % minutes of availability of cars
%         E.electricityprice=melep(1:mthor)/1000; % convert to [$/kWh]
%         E.emissionsGridProfile=mco2(1:mthor); % [g/kWh]
%         
%         maxc=E.dkav*E.maxchargeminute/aggregateratio; % max exchangeable energy per time step [kWh]
% 
%         ELayerResults=aevopti11(E);
% 
%         contracted=FCRbid(ELayerResults,maxc)
        
    else
        contracted=P.FCR.contracted;
    end 
    
    af=contracted*1000/P.Tech.battery/60*P.e;    % FCR rate per time step (normalized)
    LimitFCR=ceil(contracted*1000/P.Tech.chargekw);
else
    FCR=false;
    LimitFCR=0;
end


%% setup pricing module 

if isfield(P,'Pricing') && ~isempty(P.Pricing)
    
    % TODO: add option to have alternative specific to each passenger
    % TODO: modify model to have inbound/outbound pricing at nodes (clusters)
    
    tp=P.Pricing.tp;
    dynamicpricing=P.Pricing.dynamic;
    ParPricing.gamma_r=P.Pricing.relocationcost; % relocation cost per minute
    ParPricing.gamma_p=P.Pricing.basetariff; % base tariff per minute
    ParPricing.gamma_alt=P.Pricing.alternative;
    ParPricing.VOT=P.Pricing.VOT;
    ParPricing.pricingwaiting=P.Pricing.pricingwaiting;
    
else
    
    tp=1;
    dynamicpricing=0;
    ParPricing.gamma_r=0; % relocation cost per minute
    ParPricing.gamma_p=0; % base tariff per minute
    ParPricing.gamma_alt=0;
    ParPricing.VOT=0;
    ParPricing.pricingwaiting=1;

end

ParPricing.c=Trs*P.e;
ParPricing.relocation=AutoRelo;

% function to calculate probability of acceptance given a certain price for each OD pair
ProbAcc=@(p,s) exp(-p.*ParPricing.c-s)./(exp(-p.*ParPricing.c-s)+exp(-ParPricing.gamma_alt*ParPricing.c));
    
if P.modechoice==0
    
    Multiplier1=1;
    Multiplier2=1;
    prices=0;
    
else
    
    switch dynamicpricing

        case 0

            prices=ParPricing.gamma_p;

            Multiplier1=ProbAcc(ParPricing.gamma_p,0);
            Multiplier2=ProbAcc(ParPricing.gamma_p,0);

        case 1

            % initialize matrix of real prices offered
            prices=ones(nc,nc,ceil(tsim/tp)+1)*ParPricing.gamma_p;

        case 2

            prices=ones(ceil(tsim/tp)+1,nc*2)*ParPricing.gamma_p;
            ParPricing.gamma_d=bestp((1:120)',ParPricing.gamma_r,ParPricing.gamma_alt);
            PerDistanceTariff=ParPricing.gamma_d(max(1,Trs)).*ParPricing.c;
            
        case 3
            
            ParPricing.gamma_d=bestp((1:120)',ParPricing.gamma_r,ParPricing.gamma_alt);
            prices=ParPricing.gamma_d(max(1,Trs));
            Multiplier1=ProbAcc(prices,0);
            Multiplier2=ProbAcc(prices,0);

    end
    
end


%% parameters for trip assignment

Par=struct('Tr',Tr,'ad',ad,'e',P.e,'minsoc',P.Operations.minsoc,'modechoice',P.modechoice,...
    'maxwait',P.Operations.maxwait,'VOT',ParPricing.VOT,'WaitingCostToggle',ParPricing.pricingwaiting,'LimitFCR',LimitFCR);


%% variables for progress display and display initializations

S.starttime=cputime;
S.lasttime=S.starttime;
S.trlayerCPUtime=zeros(tsim,1);
S.enlayerCPUtime=zeros(etsim,1);
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
        IdleReached=(g.*(1-AtChargingStation)>=P.Operations.maxidle/P.e);
        ui(IdleReached)=chargingStations(Clusters(ui(IdleReached)));
        relodistCS=Tr(sub2ind(size(Tr),u(i,IdleReached),ui(IdleReached)));
        di(IdleReached)=relodistCS;
        s(3,IdleReached)=true;
        relodist(i)=relodist(i)+sum(relodistCS);
    end
    
    
    %% charging optimization
    
    if rem(i,P.beta/P.e)==1 % only cases with energy layer
        
        % current time
        lasttimemacro=cputime;
        
        % index of energy layer
        t=(i-1)/(P.beta/P.e)+1;
        
        switch P.enlayeralg
            
            case 'no'
                
                % charge as much as possible
                zmacro(:,t)=[1;0;1;0];
                
            case 'night'
                
                % charge as much as possible only between midnight and 5am
                if i*P.e<60*5
                    zmacro(:,t)=[1;0;1;0];
                    DayCharge=0;
                else
                    zmacro(:,t)=[0;0;1;0];%[(q(i,:)<0.3);0;1;0];
                    DayCharge=1;
                end
            
            case 'aggregate'
                
                % dynamic variables
                E.maxchargeminute=P.Tech.chargekw/60*aggregateratio;    % energy exchangeable per minute per vehicle [kWh]
                actualminsoc=min(P.Operations.minsoc+P.EnergyLayer.extrasoc,mean(q(i,:))*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of energy layer
                E.storagemin=P.Tech.battery*P.m*actualminsoc; % kWh

                dktripnow=Trips.dktrip(t:t+mthor-1);    % minutes spent traveling during this horizon
                E.einit=sum(q(i,:))*P.Tech.battery;     % total initial energy [kWh]
                E.etrip=dktripnow*P.Tech.consumption;   % energy used per step [kWh] 
                E.dkav=max(0,P.m*P.beta-dktripnow);         % minutes of availability of cars
                E.electricityprice=melep(t:t+mthor-1)/1000; % convert to [$/kWh]
                E.emissionsGridProfile=mco2(t:t+mthor-1); % [g/kWh]
                maxc=E.dkav*E.maxchargeminute; % max exchangeable energy per time step [kWh]

                ELayerResults=aevopti11(E);

                if ~isempty(ELayerResults)
                    
                    % make sure that there is no discharging when V2G is not allowed
                    if P.Operations.v2g==0
                        ELayerResults.discharging(:)=0;
                    end

                    % zmacro: [relative charging, relative discharging, max charging allowed, energy required by trips]
                    zmacro(:,t:t+mthor-1)=[ ELayerResults.charging , ...
                                            ELayerResults.discharging , ...
                                            maxc , ...
                                            E.etrip]';
                    zmacro(isnan(zmacro))=0;
                    
                else
                    
                    % charge as much as possible
                    zmacro(:,t:t+mthor-1)=[ones(size(maxc,1),1),zeros(size(maxc,1),1),maxc,E.etrip]';
                    
                end

            otherwise
            
                error('energy layer must be either ''no'' or ''aggregate'' ');

        end
        
        % record time
        S.enlayerCPUtime(t)=cputime-lasttimemacro;
        S.lasttime=cputime;
        
    end
    
    v2gallowed=q(i,:)>P.Operations.v2gminsoc;
    LowSocRefill=(q(i,:)<0.6)*DayCharge;
    chargevector=(ones(1,P.m)*(zmacro(1,t)/zmacro(3,t))-v2gallowed*(zmacro(2,t)/zmacro(3,t))+LowSocRefill)*ac;
    
                
    %% pricing optimization

    % current pricing number
    kp=ceil(i/tp);

    if P.modechoice
        
        if dynamicpricing>0 && dynamicpricing<3 && (i==1 || mod(i-(ts+tr+1),tp)==0)
%         if dynamicpricing>0 && mod(i-1,tp)==0 %(i==1 || mod(i-(ts+tr+1),tp)==0)

            % current pricing calculation
            PricingStep=ceil((i-1)/tp)+1;

            % expected trips
            StartTime=(PricingStep-1)*tp+1;
            Selection0=AbuckC(StartTime)+1:AbuckC(min(length(AbuckC),StartTime+tp-1));
            a_tp=sparse(As(Selection0,1),As(Selection0,2),1,nc,nc);%+q_t;

            % number of vehicles at each station (including vehicles directed there)
            ParPricing.v=histc(Clusters(ui),1:nc);
            a_tp(1:nc+1:end)=0;
            ParPricing.a=a_tp;

            if dynamicpricing==1
            
                [pricesNow,~,~]=NLPricing4(ParPricing); % OD-pair-based pricing
                
                prices(:,:,PricingStep)=pricesNow;
            
                Multiplier1=(ProbAcc(prices(:,:,kp),0));
                Multiplier2=(ProbAcc(prices(:,:,kp+1),0));
                
            elseif dynamicpricing==2
                
                [pricesNow,~,~]=NLPricing2(ParPricing); % node-based pricing
                
                prices(PricingStep,:)=pricesNow';
                % prices(PricingStep,:)=[pricesNow(1:nc)-pricesNow(nc+1:nc*2) ; pricesNow(nc*2+1:nc*3)-pricesNow(nc*3+1:nc*4)]';
                
                Surcharges1=prices(kp,1:nc)+prices(kp,nc+1:2*nc)';
                Surcharges2=prices(kp+1,1:nc)+prices(kp+1,nc+1:2*nc)';
                
                Multiplier1=(ProbAcc(PerDistanceTariff,Surcharges1));
                Multiplier2=(ProbAcc(PerDistanceTariff,Surcharges2));
                
            end

            % plot(0:0.01:0.5,histc(normalizedprices(:)*2*m.gamma_p,0:0.01:0.5))
            
        end
        
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

        % expected trips
        NextPricing=min(kp*tp+1,length(AbuckC));

        Selection1a=AbuckC(i)+1:AbuckC(min(length(AbuckC),min(NextPricing-1,i+ts)));
        Selection1b=AbuckC(NextPricing)+1:AbuckC(min(length(AbuckC),i+ts));
        a_ts=(Multiplier1.*sparse(As(Selection1a,1),As(Selection1a,2),1,nc,nc))+...
             (Multiplier2.*sparse(As(Selection1b,1),As(Selection1b,2),1,nc,nc));

        Selection2a=AbuckC(i)+1:AbuckC(min(length(AbuckC),min(NextPricing-1,i+ts+tr)));
        Selection2b=AbuckC(NextPricing)+1:AbuckC(min(length(AbuckC),i+ts+tr));
        a_to=(Multiplier1.*sparse(As(Selection2a,1),As(Selection2a,2),1,nc,nc))+...
             (Multiplier2.*sparse(As(Selection2b,1),As(Selection2b,2),1,nc,nc));
        
        % Vin: vehicles information in the form: [station delay soc connected relocating]
        Vin=[Clusters(ui) , di' , q(i,:)' , s(2,:)' , logical(s(1,:)+s(3,:))' ];
        ParRel.ad=ad;
        ParRel.dw=dw; % number of passengers waiting at each station
        ParRel.a_ts=round(sum(a_ts))'; % expected arrivals between now and now+ts
        ParRel.a_to=round(sum(a_to,2)); % expected requests between now and now+ts+tr
        ParRel.Trs=Trs;
        ParRel.bmin=bmin;
        ParRel.LimitFCR=LimitFCR;
        
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
    
        % TODO: fix the case for mode choice without pricing, harmonize code
        % TODO: fix pooling option 
        
        % calculate pricing
        switch dynamicpricing
            case 0
                pp=TripDistances(trips)*prices;
            case 1
                SelectorClusters=sub2ind(size(Trs),As(trips,1),As(trips,2));
                pricesNow=prices(:,:,kp);
                pp=pricesNow(SelectorClusters).*TripDistances(trips);
            case 2
                Surcharge=prices(kp,As(trips,1))+prices(kp,nc+As(trips,2));
                pp=(   ParPricing.gamma_d(TripDistances(trips)).*(TripDistances(trips)')  +  Surcharge  )';
            case 3
                pp=ParPricing.gamma_d(TripDistances(trips))'.*TripDistances(trips);
        end
        alte=exp(-TripDistances(trips)*ParPricing.gamma_alt);
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
        
        % power exchanged for vehicles charging
        acv=(q(i,:)<P.FCR.fastchargesoc)*ac+(q(i,:)>=P.FCR.fastchargesoc)*ac*P.FCR.slowchargeratio;
        
        %% create set points at beginning of charging period
        
        if rem(i,P.beta/P.e)==1
            % aggregate set point (kWh)
            SetPointUpPeriod=(zmacro(1,t));  % set point of aggregate fleet (kWh)
            SetPointDownPeriod=(zmacro(2,t));  % set point of aggregate fleet (kWh)
            % expected total vehicle capacity in the period (kWh)
            CapUpPeriod=max(0,s(2,:).*min(acv*P.beta/P.e,P.Operations.maxsoc-q(i,:))*P.Tech.battery); % charge
            CapDownPeriod=max(0,s(2,:).*min(acv*P.beta/P.e,(q(i,:)-P.Operations.v2gminsoc)*Efficiency)*P.Tech.battery); % discharge
            % set point for each vehicle for each time step (kWh)
            SetPointUp=min(SetPointUpPeriod,sum(CapUpPeriod))/P.beta*P.e;  % set point of aggregate fleet (kWh)
            SetPointDown=min(SetPointDownPeriod,sum(CapDownPeriod))/P.beta*P.e;  % set point of aggregate fleet (kWh)
        end
        
        
        %% charging 
        
        % available power from fleet
        maxsoceff=1; % maxsoceff=P.Operations.maxsoc;
        v2gallowed=q(i,:)>P.Operations.v2gminsoc;
        % actual capacity in this time step
        CapUp=s(2,:).*min(acv,maxsoceff-q(i,:)); % charge
        CapDown=s(2,:).*v2gallowed.*min(acv,(q(i,:)-P.Operations.minsoc)*Efficiency); % discharge
        
        % calculate ratios
        if sum(CapUp)>0
            eRatioUp=min(1,SetPointUp/(sum(CapUp)*P.Tech.battery));
        else
            eRatioUp=0;
        end
        if sum(CapDown)>0
            eRatioDown=min(1,SetPointDown/(sum(CapDown)*P.Tech.battery));
        else
            eRatioDown=0;
        end
        
        % calculate charging for each vehicle
        e(i,:)=CapUp*eRatioUp-CapDown*eRatioDown;
        
        
        %% FCR provision
        
        % needed FCR power
        FCRNeed=(f(i)-50)/(P.FCR.limits(2)-P.FCR.limits(1))*2;
        FCRNeedUp=af*min(1,max(0,FCRNeed)); % charge
        FCRNeedDown=af*min(1,max(0,-FCRNeed)); % discharge
        
        % available power from fleet
        AvailableUp=s(2,:).*min(ac-max(0,e(i,:)),1-(q(i,:)+e(i,:))); % charge
        AvailableDown=s(2,:).*min(ac-max(0,-e(i,:)),(q(i,:)+e(i,:))-P.Operations.minsoc); % discharge
        
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
        ef(i,:)=AvailableUp.*FCRRatioUp-AvailableDown.*FCRRatioDown;
           
    else
        
        % capacity based
        CapUp=s(2,:).*min(ac,P.Operations.maxsoc-q(i,:)); % charge
        CapDown=s(2,:).*min(ac,(q(i,:)-P.Operations.minsoc)*Efficiency); % discharge
        
        e(i,:)=min(CapUp,max(0,chargevector))+max(-CapDown,min(0,chargevector));
%         e(i,:)=s(2,:).*max(-ac,min(ac,(min(P.Operations.maxsoc,max( P.Operations.minsoc , q(i,:)+chargevector ) )-q(i,:))));
        
    end

    % update SOC 
    %     q(i+1,:)=min(P.Operations.maxsoc,max(P.Operations.minsoc,q(i,:)+e(i,:)-(d(i,:)>0).*ad));
    q(i+1,:)=q(i,:)+max(0,e(i,:)+ef(i,:))+min(0,e(i,:)+ef(i,:))/Efficiency-(di>0).*ad;

    % update position
    u(i+1,:)=ui;%(i,:);
    
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
    S.trlayerCPUtime(i)=cputime-S.lasttime;
    S.lasttime=cputime;
    
    % update cputime of this step
	S.comptime(i+1)=cputime;
    
end


%% final results

% main Sim struct
% if dispiter<0
%     % if 
%     uChanges=u.*logical(u(1:end,:)~=[zeros(1,size(u,2));u(1:end-1,:)]);
%     [I,J,V]=find(uChanges);
%     Sim.uSummary=uint16([I,J,V]);
% else
Sim.u=uint16(u); % final destination of vehicles (station) [tsim x m]
% end
Sim.q=single(q); % state of charge 
Sim.e=sparse(e/ac*P.Tech.chargekw);
if FCR
    Sim.ef=single(ef*P.Tech.battery*60/P.e);
else
    Sim.ef=sparse(ef);
end
Sim.waiting=sparse(waiting); % waiting times
Sim.dropped=sparse(dropped); % dropped requests
Sim.chosenmode=chosenmode; % chosen mode
Sim.waitingestimated=sparse(waitingestimated); % estimated waiting time (only mode choice)
Sim.relodist=relodist*P.e; % relocation minutes
Sim.tripdist=tripdist*P.e; % trip minutes
Sim.emissions=(sum(Sim.e/60*P.e,2)')*co2(1:tsim)/10^6; % emissions [ton]

% add pricing info
if isfield(P,'Pricing')
    
    Sim.revenues=sum((offeredprices-ParPricing.gamma_r.*TripDistances).*chosenmode.*(1-dropped));
    Sim.relocationcosts=sum(relodist)*P.e*ParPricing.gamma_r;
    Sim.offeredprices=offeredprices;
    Sim.prices=prices;
    
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
Res.cost=(sum(Sim.e/60/1000*P.e,2)')*elep(1:tsim)+Sim.emissions*P.carbonprice;
Res.dropped=sum(dropped)/length(A);
Res.peakwait=max(waiting);
Res.avgwait=mean(waiting);

% save results
if extsave>0
    save(simname,'Res','P');
end


%% end display

meanqnow=mean(q(i,:));
if dispiter<0
    fprintf('sim #%d ',-dispiter)
end
fprintf('successfully completed - avg soc: %0.2f - total time: %d:%0.2d - \n',meanqnow,floor(elapsed/60),round(rem(elapsed,60)));

end




