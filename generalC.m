% [Res]=GENERALC(P[,extsave,dispiter])
% Run SAEV simulation and relocation/charge oprimization.
% Vehicles start at beginning of time step, arrive at end of time step.
% 
% TODO: explicitly define and document all statuses and corresponding
% physical meanings
% TODO: fix problem with pickup from same node (Tr(i,i)=0?)
% 
% see also CPAR

function [Res]=generalC(P,extsave,dispiter)

if strcmp(P.trlayeralg,'opti')
    Res=generalOpti(P,extsave,dispiter);
    return
end

%% initializations

addpath functions utilities
DataFolder=setDataFolder();


%% parallel computing and input check

parcomp=is_in_parallel();
if nargin<3
    dispiter=1;
    if nargin<2
        extsave=1;
    end
end

if isfield(P,'tripfolder')
    TripName=P.tripfolder;
else
    TripName=P.tripfile;
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


%% load external files: scenario, trips and energy

% load distance matrix
load([DataFolder 'scenarios/' P.scenario '.mat'],'T');

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
ts=round(P.TransportLayer.ts/P.e);
tr=round(P.TransportLayer.tr/P.e);
tx=round(P.TransportLayer.tx/P.e);
bmin=P.TransportLayer.bmin;
mthor=round(P.EnergyLayer.mthor/P.beta);

% main variables
q=zeros(tsim,P.m,'double'); % SOC
e=zeros(tsim,P.m,'double'); % charging
u=zeros(tsim,P.m,'double'); % vehicles in charging stations
d=zeros(tsim,P.m,'double'); % delay
g=zeros(tsim,P.m,'double'); % idle time
s1=false(tsim,P.m); % relocating status
s2=false(tsim,P.m); % connection status
s3=false(tsim,P.m); % moving to charging station status
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

% parameters for trip assignment
Par=struct('Tr',Tr,'ad',ad,'e',P.e,'minsoc',P.Operations.minsoc,'modechoice',logical(P.modechoice+isfield(P,'pricing')),'maxwait',P.Operations.maxwait);

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


%% setup clustering

if isfield(P,'clusters')
    chargingStations=P.chargingStations;
    Clusters=P.clusters;
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
b=zeros(ceil(tsim/tx),nc,'double');             % imbalance


%% initial states

q(1,:)=P.Operations.initialsoc.*ones(1,P.m);      % initial state of charge
if isfield(P.Operations,'uinit')
    u(1,:)=P.Operations.uinit;
else
    u(1,:)=chargingStations(randi(nc,1,P.m));                 % initial position of vehicles
end


%% setup charging module

if strcmp(P.enlayeralg,'aggregate') 
    
    % generate aggregate trip statistics
    EMDFileName=[TripName '-' num2str(P.tripday)];
    [Trips,~,~]=generateEMD(A,Atimes,T,etsim,EMDFileName);

    % append values for next day
    if isfield(P,'tripfolder')
        P2=P;
        P2.tripday=P.tripday+1;
        [A2,Atimes2,~,~]=loadTrips(P2);
        EMDFileName=[TripName '-' num2str(P2.tripday)];
        [Trips2,~,~]=generateEMD(A2,Atimes2,T,etsim,EMDFileName);
        Trips.dktrip=[Trips.dktrip(1:48,:) ; Trips2.dktrip(1:48,:)];
    end

    % energy layer variable: static values
    E.v2g=P.Operations.v2g; % use V2G?
    E.eta=1;                % 
    E.selling=1;            % can sell to the grid?
    E.minfinalsoc=0.9;      % final SOC. This only works for optimization horizon of ~24h
    E.T=mthor;% number of time steps in energy layer
    E.cyclingcost=P.Tech.cyclingcost;                       % battery cycling cost [$/kWh]
    E.storagemax=P.Tech.battery*P.m*P.Operations.maxsoc;    % max total energy in batteries [kWh]
    E.maxchargeminute=P.Tech.chargekw/60;                   % energy exchangeable per minute per vehicle [kWh]
    E.carbonprice=P.carbonprice;                            % carbon price [$ per kg]

    % matrix of optimal control variables for energy layer
    zmacro=zeros(4,etsim+P.EnergyLayer.mthor); 
    
else 
    
    Trips=[];
    zmacro=zeros(4,etsim);
    
end


%% setup mode choice

% TODO: harmonize and generalize
VOT=15; % value of time
beta=1.5; % tortuosity of walking
% CostMinute=P.TransportLayer.basetariff;

% calculate benefit of alternative trip
WalkingSpeed=4; % km/h
if ~isempty(Distances)
    UtilityWalking=-Distances.RawDistance*beta/WalkingSpeed*VOT;
else
    UtilityWalking=0;
end
% PTSpeed=20;
% UtilityPT=-Distances.RawDistance/PTSpeed*VOT;


%% setup pricing module 

if isfield(P.TransportLayer,'tp')
    tp=P.TransportLayer.tp;
else 
    tp=1;
end

if isfield(P,'pricing')
    
    m.c=Tr*P.e;
    m.gamma_r=P.TransportLayer.relocationcost; % relocation cost per minute
    m.gamma_p=P.TransportLayer.basetariff; % base tariff per minute
    if isfield(P.TransportLayer,'alternative')
        m.gamma_alt=P.TransportLayer.alternative;
    else
        m.gamma_alt=m.gamma_p;
    end
    
    % real prices offered
    prices=ones(n,n,ceil(tsim/tp)+1)*m.gamma_p;
    
    % calculates probability of acceptance given a certain price for each OD pair
    g=@(p) exp(-p.*m.c)./(exp(-p.*m.c)+exp(-m.gamma_alt*m.c));
    
    if isfield(P,'pricingwaiting')
        WaitingCostToggle=1;
    else
        WaitingCostToggle=0;
    end
    
else
    
	prices=zeros(1,1,ceil(tsim/tp)+1);
    
    Multiplier1=1;
    Multiplier2=1;
    
end


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
    
    
    %% charging optimization
    
    if rem(i,P.beta)==1 % only cases with energy layer
        
        % current time
        lasttimemacro=cputime;
        
        % index of energy layer
        t=(i-1)/P.beta+1;
        
        switch P.enlayeralg
            
            case 'no'
                
                % charge as much as possible
                zmacro(1,t)=1;
            
            case 'aggregate'
                
                % dynamic variables
                actualminsoc=min(P.Operations.minsoc+P.EnergyLayer.extrasoc,mean(q(i,:))*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of energy layer
                E.storagemin=P.Tech.battery*P.m*actualminsoc; % kWh

                dktripnow=Trips.dktrip(t:t+mthor-1); % time steps spent traveling during this horizon
                E.einit=sum(q(i,:))*P.Tech.battery;            % total initial energy [kWh]
                E.etrip=dktripnow*P.Tech.consumption;        % energy used per step [kWh] 
                E.dkav=max(0,P.m*P.beta*P.e-dktripnow); % minutes of availability of cars
                E.electricityprice=melep(t:t+mthor-1)/1000; % convert to [$/kWh]
                E.emissionsGridProfile=mco2(t:t+mthor-1); % [g/kWh]

                ELayerResults=aevopti11(E);

                % make sure that there is no discharging when V2G is not allowed
                if P.Operations.v2g==0
                    ELayerResults.discharging(:)=0;
                end

                % zmacro: [relative charging, relative discharging, max charging allowed, energy required by trips]
                maxc=E.dkav*E.maxchargeminute; % max exchangeable energy per time step [kWh]
                zmacro(:,t:t+mthor-1)=[ ELayerResults.charging./maxc , ...
                                                                        ELayerResults.discharging./maxc , ...
                                                                        maxc , ...
                                                                        E.etrip]';
                zmacro(isnan(zmacro))=0;

            otherwise
            
                error('energy layer must be either ''no'' or ''aggregate'' ');

        end
        
        % record time
        S.enlayerCPUtime(t)=cputime-lasttimemacro;
        S.lasttime=cputime;
        
    end
    
    v2gallowed=q(i,:)>P.Operations.v2gminsoc;
    chargevector=(ones(1,P.m)*zmacro(1,t)-v2gallowed*zmacro(2,t))*ac;
    
                
    %% pricing

    % current pricing number
    kp=ceil(i/tp);

    if isfield(P,'pricing')  % TODO: waiting trips are not influenced by prices! should be included outside

        if P.pricing~=0 && (i==1 || mod(i-(ts+tr+1),tp)==0)

            % number of vehicles at each station (including vehicles directed there)
            % uv=histc(u(i,:)+sum(v(i:end,:))+sum(w(i:end,:)),1:n);
            uv=histc(u(i,:),1:n);

            % current pricing calculation
            PricingStep=ceil((i-1)/tp)+1;

            % expected trips
            StartTime=(PricingStep-1)*tp+1;
            Selection0=AbuckC(StartTime)+1:AbuckC(min(length(AbuckC),StartTime+tp-1));
            a_tp=sparse(A(Selection0,1),A(Selection0,2),1,n,n);%+q_t;

            m.v=uv';
            m.a=a_tp;

            [pricesNow,~,~]=NLPricing(m);

            prices(:,:,PricingStep)=pricesNow;

            % plot(0:0.01:0.5,histc(normalizedprices(:)*2*m.gamma_p,0:0.01:0.5))

        end

        Multiplier1=(g(prices(:,:,kp)));
        Multiplier2=(g(prices(:,:,kp+1)));
    end
    

    %% relocation

    % if it's time for a relocation decision
    if mod(i-1,tx)==0
        
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
        Vin=[Clusters(u(i,:)) , d(i,:)' , q(i,:)' , s2(i,:)' , logical(s1(i,:)+s3(i,:))' ];
        Par.dw=dw; % number of passengers waiting at each station
        Par.a_ts=round(sum(a_ts))'; % expected arrivals between now and now+ts
        Par.a_to=round(sum(a_to,2)); % expected requests between now and now+ts+tr
        Par.Trs=Trs;
        Par.bmin=bmin;
        
        [Vout,bkt]=Relocation(Vin,Par);
        
        % update vehicles position
        used=logical(Vout(:,2));
        relodestinations=chargingStations(Vout(used,1));
        relodisti=Tr(sub2ind(size(Tr),u(i,used),relodestinations'));
        
        u(i,used)=relodestinations; 
        d(i,used)=d(i,used)+relodisti;
        
        % update vehicles status (relocating vehicles cannot be relocated)
        s1(i,used)=1;
        s2(i,used)=0;
        s3(i,used)=0;
        
        % update results
        b(kt,:)=bkt;
        relodist(i)=sum(relodisti);
        
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
        Selector=sub2ind(size(Tr),A(trips,1),A(trips,2));
        if isfield(P,'pricing')
            pp=prices(Selector,kp);
            alte=exp(-Tr(Selector)*P.e*m.gamma_alt);
        else
            pp=ones(length(trips),1)*prices(1,1,kp);
            alte=zeros(length(trips),1);
            %             Cost=(VOT/60+CostMinute);
            %             UtilityAlternative=exp(UtilityWalking(tripID));
        end
        offeredprices(trips)=pp;
        
        % Vin: vehicles information in the form: [station delay soc connected]
        Vin=[u(i,:)' , d(i,:)' , q(i,:)' , s2(i,:)' ];
        
        % Bin: passengers info in the form: [O D waiting offeredprice utilityalternative]
        Bin=[A(trips,:) , waiting(trips) , pp , alte ];

        % tripAssignment (no clustering) or tripAssignment2 (clustering)
        if isfield(P,'clusters')
            [Vout,Bout,tripdisti,relodistiPU,queuei]=tripAssignment2(Vin,Bin,Par);
        else 
            [Vout,Bout,tripdisti,relodistiPU,queuei]=tripAssignment(Vin,Bin,Par);
        end
        
        % update vehicles positions
        u(i,:)=Vout(:,1)';
        d(i,:)=Vout(:,2)';
        
        % update vehicles status
        used=logical(Vout(:,3));
        s1(i,used)=0;
        s2(i,used)=0;
        s3(i,used)=0;
        
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
    
    % power exchanged for vehicles charging
    e(i,:)=min(P.Operations.maxsoc,max(P.Operations.minsoc,q(i,:)+s2(i,:).*chargevector))-q(i,:);

    % update SOC 
    %     q(i+1,:)=min(P.Operations.maxsoc,max(P.Operations.minsoc,q(i,:)+e(i,:)-(d(i,:)>0).*ad));
    q(i+1,:)=q(i,:)+e(i,:)-(d(i,:)>0).*ad;

    % update position
    u(i+1,:)=u(i,:);
    
    % update delay
    d(i+1,:)=max(0,d(i,:)-1);
    
    % update idle time
    g(i+1,:)=(d(i,:)==0).*(g(i,:)+1);
    
    AtChargingStation=sum(u(i+1,:)==chargingStations);
    
    % update statuses
    s1(i+1,:)=(d(i+1,:)>0).*s1(i,:);
    s2(i+1,:)=AtChargingStation.*(d(i+1,:)==0);
    s3(i+1,:)=(d(i+1,:)>0).*s3(i,:);
    
    if isfield(P,'clusters')
        % move idle vehicles back to charging stations
        IdleReached=(g(i+1,:).*(1-AtChargingStation)>=P.Operations.maxidle/P.e);
        u(i+1,IdleReached)=chargingStations(Clusters(u(i+1,IdleReached)));
        relodistCS=Tr(sub2ind(size(Tr),u(i,IdleReached),u(i+1,IdleReached)));
        d(i+1,IdleReached)=relodistCS;
        s3(i+1,IdleReached)=true;
        relodist(i)=relodist(i)+sum(relodistCS);
    end
    
    % record time
    S.trlayerCPUtime(i)=cputime-S.lasttime;
    S.lasttime=cputime;
    
    % update cputime of this step
	S.comptime(i+1)=cputime;
    
end


%% final results

% main Sim struct
Sim.u=uint8(u); % final destination of vehicles (station) [tsim x m]
Sim.q=single(q); % state of charge 
Sim.e=sparse(e/ac*P.Tech.chargekw);
Sim.waiting=sparse(waiting); % waiting times
Sim.dropped=sparse(dropped); % dropped requests
Sim.chosenmode=chosenmode; % chosen mode
Sim.waitingestimated=sparse(waitingestimated); % estimated waiting time (only mode choice)
Sim.relodist=relodist*P.e; % relocation minutes
Sim.tripdist=tripdist*P.e; % trip minutes
Sim.emissions=(sum(Sim.e/60*P.e,2)')*co2(1:tsim)/10^6; % emissions [ton]

% add pricing info
if isfield(P,'pricing')
    
    distances=Tr(sub2ind(size(Tr),A(:,1),A(:,2)))*P.e; % minutes
    
    Sim.revenues=sum((offeredprices-m.gamma_r).*distances.*chosenmode.*(1-dropped));
    Sim.relocationcosts=sum(relodist)*P.e*m.gamma_r;
    Sim.offeredprices=offeredprices;
    Sim.prices=prices;
    
end

% Internals struct
Internals.b=b;
Internals.d=sparse(d);
Internals.g=sparse(g);
Internals.s1=sparse(logical(s1));
Internals.s2=sparse(logical(s2));
Internals.s3=sparse(logical(s3));
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




