% [Res]=GENERALC(P[,extsave,dispiter])
% Run SAEV simulation and relocation/charge oprimization with coordinates.
% Vehicles start at beginning of time step, arrive at end of time step.
%
% Mode choice - 0: walking, 1: SAEV
% 
% see also CPAR

function [Res]=generalC(P,extsave,dispiter)


%% initializations

addpath functions
addpath utilities
DataFolder=setDataFolder();


%% parallel computing and input check

parcomp=is_in_parallel();
if nargin<3
    dispiter=1;
    if nargin<2
        extsave=1;
    end
end

if strcmp(P.enlayeralg,'opti') && strcmp(P.trlayeralg,'simplified')
    warning('impossible combination!');
    return
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

load(['data/scenarios/' P.scenario '.mat'],'T','C');

% NOTE: can add secondary trip file (real vs expected/forecasted)
[A,Atimes,AbuckC,Distances]=generateGPStrips(P);
AbuckC=AbuckC(1:P.e:end);

% NOTE: should generalize vector length for cases with different beta, e,
% etc. Also: change names of variables
% elep is in $/MWh
load(['data/eleprices/' P.gridfile '.mat'],'x','y');
melep=repelem([x(:,P.gridday);x(:,rem(P.gridday,size(x,2))+1)],2/P.e,1); % electricity price profiles
if exist('y','var') % carbon emissions profiles [g/kWh]
    mco2=repelem([y(:,P.gridday);y(:,rem(P.gridday,size(y,2))+1)],2/P.e,1);
else
    mco2=zeros(24*2*2,1);
end
clear x y;


%% parameters of simulation

% parameters
n=size(T,1);              % number of nodes
tsim=1440/P.e;            % number of time steps in transport layer
etsim=tsim/P.beta;        % number of time steps in energy layer
Tr=max(1,round(T/P.e));   % distance matrix in transport layer steps
ac=round(P.Tech.chargekw/P.Tech.battery/60*P.e,3);    % charge rate per time step (normalized)
ad=P.Tech.consumption/P.Tech.battery*P.e;             % discharge rate per time step (normalized)
elep=repelem(melep,P.beta);                 % electricity price in each transport layer time step
co2=repelem(mco2,P.beta);                 % carbon intensity in each transport layer time step

% main variables
q=zeros(tsim,P.m,'double');            % SOC
e=zeros(tsim,P.m,'double');            % charging
u=zeros(tsim,P.m,'double');            % vehicles in charging stations
v=zeros(tsim+100,P.m,'double');        % auxiliary variable to assign vehicles to future for trips with passengers
w=zeros(tsim+100,P.m,'double');        % auxiliary variable to assign vehicles to future for relocation
b=zeros(etsim,n,'double');             % imbalance
s=zeros(tsim+100,n,'double');          % auxiliary variable to keep track of vehicles arriving at stations

% working variables
queue=zeros(100,1);          % temporary variable to store queued arrivals
% poolid=0;                    % index of pool

% results variables
waiting=zeros(length(A),1);  % minutes waited for each request
dropped=zeros(length(A),1);  % request is dropped?
chosenmode=false(length(A),1);% which mode is chosen?
pooling=zeros(length(A),1);  % pool ID of each user (if ride shared)
% traveled=zeros(length(A),1); % trip length (minutes)
relodist=zeros(ceil(tsim),1); % distances of relocation (at moment of decision)
tripdist=zeros(ceil(tsim),1); % distances of trips (at moment of acceptance)
waitingestimated=zeros(length(A),1);  % estimated minutes to wait for each request
offeredprices=ones(length(A),1);  % price offered to each passenger

% initial states
q(1,:)=P.Operations.initialsoc.*ones(1,P.m);      % initial state of charge
if isfield(P.Operations,'uinit')
    u(1,:)=P.Operations.uinit;
else
    u(1,:)=randi(n,1,P.m);                 % initial position of vehicles
end

%% trip processing

% if isfield(P,'pricing') && P.pricing==false
% 
%     [A,I]=tripAcceptance(A,0.5);
%     Atimes=Atimes(I,:);
%     
% end

% generate number of arrivals at each station
% generate EMD in case of aggregate energy layer
% calculate from expected arrivals
% etsim is the number of energy layer time steps in a day
% dkemd, dkod, dktrip are the number of minutes of travel for
% relocation, serving trips, and total, respectively, for each
% energy layer time step. fk
[fo,fd,Trips]=generateEMD(A,Atimes,T,etsim,TripName,P.tripday);

% calculate benefit of alternative trip
VOT=15; % value of time
MaxHor=min(15,ceil(P.Operations.maxwait/P.e)); % maximum horizon for waiting time estimation
CostMinute=0.2;
WalkingSpeed=4;
UtilityWalking=-Distances.RawDistance*1.5/WalkingSpeed*VOT;
% PTSpeed=20;
% UtilityPT=-Distances.RawDistance/PTSpeed*VOT;

%% setup energy layer

if strcmp(P.enlayeralg,'aggregate') 

    % append values for next day
    if isfield(P,'tripfolder') 
        P2=P;
        P2.tripday=P.tripday+1;
        [A2,Atimes2,~,~]=generateGPStrips(P2);
        [fo2,fd2,Trips2]=generateEMD(A2,Atimes2,T,etsim,TripName,P.tripday+1);
        fo=[fo(1:1440,:) ; fo2(1:1440,:)];
        fd=[fd(1:1440,:) ; fd2(1:1440,:)];
        Trips.dktrip=[Trips.dktrip(1:48,:) ; Trips2.dktrip(1:48,:)];
    end

    % energy layer variable: static values
    E.v2g=P.Operations.v2g; % use V2G?
    E.eta=1;                % 
    E.selling=1;            % can sell to the grid?
    E.minfinalsoc=0.9;      % final SOC. This only works for optimization horizon of ~24h
    E.T=P.EnergyLayer.mthor;% number of time steps in energy layer
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


%% setup opti transport layer

if strcmp(P.trlayeralg,'opti') 
    
    if isfinite(P.Operations.maxwait) && P.Operations.maxwait>=0
        warning('maxwait is ignored for opti transport layer: no dropped requests allowed')
    end
    
    maxt=max(Tr(:));
    varno=n^2+P.m*n*(maxt+2)+P.m;   % passengers waiting, cars positions, cars waiting, SOC
    ctrno=n^2*P.m*2+P.m*2;          % vehicles movements with passenger, for relocaiton, charging power

    % create transport layer matrices
    [namesim]=generatematrices3(n,P.m,Tr,maxt,ac,ad,P.TransportLayer.thor,P.Operations.maxsoc,P.Operations.minsoc);
    load(namesim,'Aopti','Bopti','Adis','bdis','bdisc','bdisC','Aeq','beq','beqc','intcon','lb','ub','fu','fx','fq','fqv2g','fsoc');
    
    % initializations of state variables
    uinit=zeros(n,P.m);           % vehicles waiting at a station - binary variable
    uinit(u(1,:)+(0:P.m-1)*n)=1;
    uinit=reshape(uinit,[n*P.m,1]);
    x=[zeros(n^2,1);zeros(n*P.m*(maxt+1),1);uinit;q(1,:)'];
    
    % initializations of simulation variables
    cname=['data/temp/c-' TripName '-' num2str(P.tripday) '.mat'];
    if exist(cname,'file')
        load(cname,'c');
    else
        c=double(convertAtimes(A,Atimes,n));
        save(cname,'c');
    end
    c1=cat(3,c,zeros(n,n,tsim)); % add padding
    
    % sum arrivals (counted by minute) over one time step
    if P.e>1
        L=floor(size(c1,3)/2);
        c1=permute(squeeze(sum(reshape(permute(c1(:,:,1:L*P.e),[3,1,2]),[P.e,L,n,n]),1)),[2,3,1]);
        % NOTE: implementation for stochastic arrivals needed:
        % c2=permute(squeeze(sum(reshape(permute(c2(:,:,1:L*P.e),[3,1,2]),[P.e,L,n,n]),1)),[2,3,1]);
    end
    
    % create arrival vectors for optimization (cexpected) and simulation (c)
    c=reshape(c1(:,:,1:tsim+P.TransportLayer.thor),[n^2*(tsim+P.TransportLayer.thor),1]);
    if P.mpcpredict==true
        cexpected=c;
    else
        % NOTE: need to implement case with imperfect predictions
        error('mpcpredict==false not implemented'); 
    end
    
    % create optimization variables
	X=[x zeros(length(x),tsim)];% matrix of results
    z=zeros(ctrno,tsim);        % matrix of optimal control variables
    B=zeros(varno,tsim);        % matrix of static values
    B(1:n^2,:)=reshape(c(1:n^2*tsim),n^2,tsim); % add arrivals
    
    % intlinprog options
    if parcomp || dispiter==0
        displayopt='none';
    else
        displayopt='iter';
    end
    
    % Matlab version check
    try
        options=optimoptions('intlinprog','RelativeGapTolerance',0.02,'IntegerTolerance',1e-4,'MaxTime',1000,'LPMaxIterations',5e5,'Display',displayopt);%,'MaxNodes',10000,'Heuristics','none');%,'IntegerTolerance',1e-3);
    catch
        options=optimoptions('intlinprog','TolGapRel',0.02,'MaxTime',1000,'LPMaxIter',5e5,'Display',displayopt); % matlab 2015
    end
    
    % objective function
    f1=fx+P.TransportLayer.rho1*fu ; 
    f=(f1-(fq-fqv2g)*P.TransportLayer.rho4)/P.TransportLayer.thor;
    
    % ub selector for charging constraints
    ubChargingSelector=logical(repmat( [zeros(n*n*P.m*2,1)  ;  ones(P.m*2,1) ] , P.TransportLayer.thor,1));
    
    % how many steps of energy layer are considered in transport layer (MPC)
    EMaxHorizon=ceil(P.TransportLayer.thor/P.beta);
    % EMaxHorizon=P.EnergyLayer.mthor;
    
end


%% pricing 

if isfield(P,'pricing') 
    
    m.c=Tr*P.e;
    m.gamma_r=P.TransportLayer.relocationcost; % relocation cost per minute
    m.gamma_p=P.TransportLayer.basetariff; % base tariff per minute
    
    % real prices offered
%     normalizedprices=(m.gamma_r+m.gamma_p*2)/(4*m.gamma_p);
    prices=ones(n,n,ceil(tsim/P.TransportLayer.tp)+1)*m.gamma_p;
    
else
    
	prices=zeros(n,n,ceil(tsim/P.TransportLayer.tp)+1);
    
end

% calculates probability of acceptance given a certain price for each OD pair
g=@(s) exp(-s.*m.c)./(exp(-s.*m.c)+exp(-m.gamma_p*m.c));


%% variables for progress display and display initializations

S.starttime=cputime;
S.lasttime=S.starttime;
S.trlayerCPUtime=zeros(tsim,1);
S.enlayerCPUtime=zeros(etsim,1);

% set up time variables for cputime calculation
comptime=[cputime;zeros(tsim,1)];


%% start of iterations

for i=1:tsim
    
    
	% display progress    
    
    displayState(i,tsim,dispiter,comptime(i)-comptime(1),40)
    
    
    %% energy layer
    
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

                dktripnow=Trips.dktrip(t:t+P.EnergyLayer.mthor-1); % time steps spent traveling during this horizon
                E.einit=sum(q(i,:))*P.Tech.battery;            % total initial energy [kWh]
                E.etrip=dktripnow*P.Tech.consumption;        % energy used per step [kWh] 
                E.dkav=max(0,P.m*P.beta*P.e-dktripnow); % minutes of availability of cars
                E.electricityprice=melep(t:t+P.EnergyLayer.mthor-1)/1000; % convert to [$/kWh]
                E.emissionsGridProfile=mco2(t:t+P.EnergyLayer.mthor-1); % [g/kWh]

                ELayerResults=aevopti11(E);

                % make sure that there is no discharging when V2G is not allowed
                if P.Operations.v2g==0
                    ELayerResults.discharging(:)=0;
                end

                % zmacro: [relative charging, relative discharging, max charging allowed, energy required by trips]
                maxc=E.dkav*E.maxchargeminute; % max exchangeable energy per time step [kWh]
                zmacro(:,t:t+P.EnergyLayer.mthor-1)=[ ELayerResults.charging./maxc , ...
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
    
    
    %% transport layer
    
    switch P.trlayeralg
        
        case 'opti'
            
            % adjust values of known parameters (vector b): multiply by current Param.x and add static values
            bdist=bdis*X(1:varno,i)  +  bdisc  + ...
                bdisC * [  c(n^2*(i-1)+1:n^2*i)  ; cexpected(n^2*(i)+1:n^2*(i+P.TransportLayer.thor-1))]; % bdist must be positive
            beqt=beq*X(1:varno,i)    +beqc   ;%+beqC*c(n^2*(i-1)+1:n^2*(i+P.TransportLayer.thor));
                       
            % determine which vehicles can discharge to grid
            v2gallowed=[ones(P.m,1);(q(i,:)'>P.Operations.v2gminsoc)]*ones(1,EMaxHorizon);
            
            % create charge vector
            chargevector=repmat(     reshape(repelem(zmacro(1:2,t:t+EMaxHorizon-1),P.m,1).*v2gallowed,EMaxHorizon*2*P.m,1),P.beta,1);
            
            % apply charging constraints
            ub(ubChargingSelector)=chargevector(1:P.TransportLayer.thor*P.m*2)*ac;
            
            Aequ=Aeq;
            beqtu=beqt;
        
            % transport layer optimization
            zres=intlinprog(f,intcon,Adis,bdist,Aequ,beqtu,lb,ub,options);
            
            % optimal control variables in this time step (round binary values)
            z(:,i)=[round(zres(1:P.m*n^2*2));zres(P.m*n^2*2+1:ctrno)];
            
            % calculate new Param.x (using first time step of solution)
            X(1:varno,i+1)=round(Aopti*X(1:varno,i)+Bopti*z(:,i)+B(:,i),4);
            
            % calculate new SOC
            q(i+1,:)=(X(n^2+n*P.m*(maxt+2)+1:n^2+n*P.m*(maxt+2)+P.m,  i+1  ));
            
        case {'simplified' , 'no'}       % simplified relocation 

            % current pricing number
            kp=ceil(i/P.TransportLayer.tp);
                
            %% pricing
            if isfield(P,'pricing') 
                

                if P.pricing==true && (i==1 || mod(i-(P.TransportLayer.ts+P.TransportLayer.tr+1),P.TransportLayer.tp)==0)

                    % number of vehicles at each station (including vehicles directed there)
                    uv=histc(u(i,:)+sum(v(i:end,:))+sum(w(i:end,:)),1:n);
                    
                    % current pricing calculation
                    PricingStep=ceil((i-1)/P.TransportLayer.tp)+1;
                    
                    % NOTE: waiting trips are not influenced by prices! should be moved outside

                    % expected trips
                    StartTime=(PricingStep-1)*P.TransportLayer.tp+1;
                    %q_t=sparse(A(queue(queue>0),1),A(queue(queue>0),2),1,n,n);
                    Selection0=AbuckC(StartTime)+1:AbuckC(min(length(AbuckC),StartTime+P.TransportLayer.tp-1));
                    a_tp=sparse(A(Selection0,1),A(Selection0,2),1,n,n);%+q_t;

                    m.v=uv';
                    m.a=a_tp;

                    % [x0,pricesNow]=RelocationPricing2(m);
                    [pricesNow,~,~]=NLPricing(m);

                    prices(:,:,PricingStep)=pricesNow;
                    
                    % plot(0:0.01:0.5,histc(normalizedprices(:)*2*m.gamma_p,0:0.01:0.5))
                    
                end
            end
            
            
            %% relocation

            % if it's time for a relocation decision
            if mod(i-1,P.TransportLayer.tx)==0 && strcmp(P.trlayeralg,'simplified')

                % current relocation number
                kt=(i-1)/P.TransportLayer.tx+1;
                
                % number of vehicles at each station
                uv=histc(u(i,:),1:n);

                % number of waiting passenger at station
                if sum(queue)>0
                    dw=histc(A(queue(queue>0),1)',1:n);
                else 
                    dw=zeros(1,n);
                end
                
%                 % expected trips
%                 Selection1=AbuckC(i)+1:AbuckC(min(length(AbuckC),i+P.TransportLayer.ts));
%                 a_ts=round((1-prices(:,:,kp)).*sparse(A(Selection1,1),A(Selection1,2),1,n,n));
%                 Selection2=AbuckC(i)+1:AbuckC(min(length(AbuckC),i+P.TransportLayer.ts+P.TransportLayer.tr));
%                 a_to=round((1-prices(:,:,kp)).*sparse(A(Selection2,1),A(Selection2,2),1,n,n));
                
                % expected trips
                NextPricing=min(ceil(i/P.TransportLayer.tp)*P.TransportLayer.tp+1,length(AbuckC));
                
                Selection1a=AbuckC(i)+1:AbuckC(min(length(AbuckC),min(NextPricing-1,i+P.TransportLayer.ts)));
                Selection1b=AbuckC(NextPricing)+1:AbuckC(min(length(AbuckC),i+P.TransportLayer.ts));
%                 a_ts=round((1-prices(:,:,kp)).*sparse(A(Selection1a,1),A(Selection1a,2),1,n,n))+...
%                      round((1-prices(:,:,kp+1)).*sparse(A(Selection1b,1),A(Selection1b,2),1,n,n));
                a_ts=((g(prices(:,:,kp))).*sparse(A(Selection1a,1),A(Selection1a,2),1,n,n))+...
                     ((g(prices(:,:,kp+1))).*sparse(A(Selection1b,1),A(Selection1b,2),1,n,n));
                
                Selection2a=AbuckC(i)+1:AbuckC(min(length(AbuckC),min(NextPricing-1,i+P.TransportLayer.ts+P.TransportLayer.tr)));
                Selection2b=AbuckC(NextPricing)+1:AbuckC(min(length(AbuckC),i+P.TransportLayer.ts+P.TransportLayer.tr));
%                 a_to=round((1-prices(:,:,kp)).*sparse(A(Selection2a,1),A(Selection2a,2),1,n,n))+...
%                      round((1-prices(:,:,kp+1)).*sparse(A(Selection2b,1),A(Selection2b,2),1,n,n));
                a_to=((g(prices(:,:,kp))).*sparse(A(Selection2a,1),A(Selection2a,2),1,n,n))+...
                     ((g(prices(:,:,kp+1))).*sparse(A(Selection2b,1),A(Selection2b,2),1,n,n));
                
                % expected imbalance at stations
                b(kt,:)=uv ...
                    -dw ...  number of passengers waiting at each station
                    +round(sum(a_ts)) ...  expected arrivals between now and now+ts
                    -round(sum(a_to,2))' ...     expected requests between now and now+ts+tr
                    +histc(reshape(w(i:i+P.TransportLayer.ts,:),P.m*(P.TransportLayer.ts+1),1),1:n)';% vehicles relocating here between now and now+ts

%                     b(kt,:)=uv ...
%                         -dw ...  number of passengers waiting at each station
%                         +sum(fd((i-1)*P.e+1:(i+P.TransportLayer.ts)*P.e,:)) ...  expected arrivals between now and now+ts
%                         -sum(fo((i-1)*P.e+1:(i+P.TransportLayer.ts+P.TransportLayer.tr)*P.e,:)) ...     expected requests between now and now+ts+tr
%                         +histc(reshape(w(i:i+P.TransportLayer.ts,:),P.m*(P.TransportLayer.ts+1),1),1:n)';% vehicles relocating here between now and now+ts

                % identify feeder and receiver stations
                F=min(uv,(b(kt,:)-P.TransportLayer.bmin).*(b(kt,:)>=P.TransportLayer.bmin)); % feeders
                R=(-b(kt,:)+P.TransportLayer.bmin).*(b(kt,:)<P.TransportLayer.bmin); % receivers

                % identify optimal relocation flux
                x=optimalrelocationfluxes(F,R,Tr);
                
                if ~isempty(x)

                    % read results
                    [Fs,Rs,Vr]=find(x);

                    % distance of relocation
                    arri=Tr(sub2ind(size(Tr),Fs,Rs)); 

                    % duplicate fluxes with multiple vehicles
                    Fs=repelem(Fs,Vr);
                    Rs=repelem(Rs,Vr);
                    arri=repelem(arri,Vr);

                    % satisfy longer relocation tasks first
                    [arris,dstnid]=sort(arri,'descend');

                    for ka=1:length(arris)

                        % find candidate vehicles for the task with enough soc
                        candidates=q(i,:).*(u(i,:)==Fs(dstnid(ka))).*(q(i,:)/ad >= arris(ka)); 

                        % remove unavailable vehicles
                        candidates(candidates==0)=NaN;

                        % sort candidate vehicles by SOC
                        [ur,ui]=max(candidates);

                        % if there is a vehicle available
                        if ~isnan(ur)

                            % update destination station
                            u(i,ui)=0;

                            % send relocation instruction
                            w(i+arris(ka),ui)=Rs(dstnid(ka)); % should be -1? depends if I assume that it starts at beginning of time period or not. Need to be explicit

                            % keep track of vehicles arriving at stations
                            s(i+arris(ka),Rs(dstnid(ka)))=s(i+arris(ka),Rs(dstnid(ka)))+1;
                            
                            % save length of relocation
                            relodist(i)=relodist(i)+arris(ka);

                        end
                    end
                end
            end
            
            
            %% charging variables

            v2gallowed=q(i,:)>P.Operations.v2gminsoc;
            chargevector=(ones(1,P.m)*zmacro(1,t)-v2gallowed*zmacro(2,t))*ac;
            
        
            %% trip assignment

            % generate trip requests for this time step
            % trips=(AbuckC((i-1)*P.e+1)+1:AbuckC(i*P.e+1))';
            trips=(AbuckC(i)+1:AbuckC(i+1))';

            % initialize 
            ql=0;

            % add queued requests
            if sum(queue)>0
                trips=[queue(queue>0);trips];
                queue(:)=0;
            end

            % if there are trips
            if ~isempty(trips)

                % for each station
                for j=1:n

                    % trips starting at station k
                    tripsK=trips((A(trips,1)==j));

                    % if there are passengers waiting at station
                    if ~isempty(tripsK)

                        % vehicles at this station
                        uid=find(u(i,:)==j);
                        
                        if P.modechoice
                            % how many vehicles have destination this station in next 20 minutes
                            % DirectedHere=sum(v(i:i+MaxHor,:)==j,2)+sum(w(i:i+MaxHor,:)==j,2);
                            % DirectedHereSum=cumsum(DirectedHere);
                            DirectedHereSum=cumsum(s(i:i+MaxHor,j));
                        end
                        
                        % soc of vehicles at this station (sorted by high soc)
                        [qj,usortid]=sort(q(i,uid),'descend');

                        % vehicle ID:
                        % uid(usortid)

                        % destination station
                        destinations=A(tripsK,2);

                        % distance of each trip
                        distancetomove=Tr(j,destinations);
                        
                        % trip priority: highest waiting first, then longest travel time
                        trippriority=distancetomove+waiting(tripsK)'*max(Tr(:));

                        % sort trips by distance (highest first)
                        [~,sortid]=sort(trippriority,'descend');
                        
                        % sort trip distances
                        distancetomovesorted=distancetomove(sortid);

                        % for each trip
                        for ka=1:length(distancetomovesorted)
                            
                            % trip ID
                            tripID=tripsK(sortid(ka));
                            
                            % candidates
                            candidates=(qj>=distancetomovesorted(ka)*ad+P.Operations.minsoc);

                            % is there a vehicle available?
                            VehicleAvailable=(sum(candidates)>0);
                            
                            if P.modechoice
                            
                                % avoid changing chosen mode after deciding
                                if chosenmode(tripID)==0

                                    % NOTE: implement choice model here
                                    % 1. estimate waiting time
                                    % 2. estimate cost
                                    % 3. decide 

                                    if VehicleAvailable

                                        WaitingTime=0;

                                    else

                                        FirstAvailable=find(DirectedHereSum>=ka,1);

                                        % if not in u, v, w, == inf
                                        if isempty(FirstAvailable)
                                            WaitingTime=Inf;
                                        else
                                            WaitingTime=FirstAvailable*P.e;
                                        end

                                    end

                                    UtilitySAEV=-distancetomovesorted(ka)*P.e*(VOT/60+CostMinute)-WaitingTime*VOT/60;
                                    AcceptProbability=exp(UtilitySAEV)/(exp(UtilitySAEV)+exp(UtilityWalking(tripID)));

                                    chosenmode(tripID)=(rand()<AcceptProbability);

                                    waitingestimated(tripID)=WaitingTime;
                                    
                                end
                                
                            else
                                
                                % NOTE: temporary solution
                                
                                if isfield(P,'pricing')
                                
                                    % avoid changing chosen mode after deciding
                                    if chosenmode(tripID)==0
                                    
                                        offeredprices(tripID)=prices(A(tripID,1),A(tripID,2),kp);

%                                         chosenmode(tripID)=(rand()>offeredprices(tripID));
                                        chosenmode(tripID)=(rand() < ...
                                            exp(-offeredprices(tripID)*distancetomovesorted(ka)*P.e) / ...
                                          ( exp(-offeredprices(tripID)*distancetomovesorted(ka)*P.e) + exp(-distancetomovesorted(ka)*P.e*m.gamma_p)));
                                        
                                    end
                                    
                                else
                                
                                    chosenmode(tripID)=1;
                                
                                end
                               
                            end
                            
                            if chosenmode(tripID)==1

                                % if there are vehicles with enough soc
                                if VehicleAvailable

                                    % vehicle with highest soc goes
                                    usortedi=find(candidates,1);

                                    % vehicle id
                                    ui=uid(usortid(usortedi));

                                    % accept request and update vehicle position
                                    u(i,ui)=0;
                                    qj(usortedi)=0;
                                    v(i+distancetomovesorted(ka),ui)=destinations(sortid(ka));
                                    
                                    % keep track of vehicles arriving at stations
                                    s(i+distancetomovesorted(ka),destinations(sortid(ka)))=s(i+distancetomovesorted(ka),destinations(sortid(ka)))+1;
                                    
                                    % update travelled distance
                                    tripdist(i)=tripdist(i)+distancetomovesorted(ka);

                                else

                                    % check if wait time of this request is less than max. wait time (minutes)
                                    if waiting(tripID)<P.Operations.maxwait

                                        % increase waiting time (minutes) for this trip
                                        waiting(tripID)=waiting(tripID)+P.e;

                                        % add this trip to the queue
                                        ql=ql+1;  % current trip
                                        queue(ql)=tripID;

                                    else

                                        % if max waiting exceeded, request dropped
                                        if pooling(tripID)>0
                                            tripsdropped=find(pooling==pooling(tripID));
                                        else
                                            tripsdropped=tripID;
                                        end

                                        % register this as a dropped request
                                        dropped(tripsdropped)=1;

                                    end

                                end
                            else
                                
                                % walking
                                
                            end
                        end
                    end
                end
            end


            %% simulation variables update

            % update SOC for vehicles charging
            e(i,:)=min(P.Operations.maxsoc,max(P.Operations.minsoc,q(i,:)+(u(i,:)>0).*chargevector))-q(i,:);

            % update SOC 
            q(i+1,:)=min(P.Operations.maxsoc,max(P.Operations.minsoc,q(i,:)+e(i,:)-(u(i,:)==0).*ad));

            % update idle vehicle positions
            u(i+1,:)=u(i,:)+v(i,:)+w(i,:);
            
        otherwise
            
            error('transport layer must be either ''opti'' or ''simplified'' ');

    end
    
    % record time
    S.trlayerCPUtime(i)=cputime-S.lasttime;
    S.lasttime=cputime;
    
    % update cputime of this step
	comptime(i+1)=cputime;
    
end


%% final calculations

if strcmp(P.trlayeralg,'simplified') 
    
    Internals.b=b;
    Internals.v=sparse(v);
    Internals.w=sparse(w);
    
end

if strcmp(P.trlayeralg,'opti') 
    
    Internals.X=sparse(X);
    Internals.z=sparse(z);
    
    d=X(1:n^2,:);
    %p=X(n^2+1:n^2+n*P.m*(maxt+1),:);
    u0=X(n^2+n*P.m*(maxt+1)+1:n^2+n*P.m*(maxt+1)+n*P.m,:);
    
    [waiting,waitinginfo]=calcwaitingtimes(P.e,c,d);
    
    Sim.waitinginfo=waitinginfo;
    
    e=z(P.m*n^2*2+1:P.m*n^2*2+P.m,:)'-z(P.m*n^2*2+P.m+1:end,:)';
    
    convertMatrix=kron(eye(P.m),(1:n)');
    u=(u0')*convertMatrix;
    
    relodist=(z(P.m*n^2+1:P.m*n^2*2,:)')*repmat(reshape(Tr,n^2,1),P.m,1);
    
end

Internals.zmacro=zmacro;

Sim.u=uint8(u); % final destination of vehicles (station) [tsim x m]
Sim.q=single(q); % state of charge 
Sim.e=sparse(e/ac*P.Tech.chargekw);

% waiting times
Sim.waiting=sparse(waiting);

% dropped requests
Sim.dropped=sparse(dropped);

% chosen mode
Sim.chosenmode=chosenmode;

% estimated waiting time (only mode choice)
Sim.waitingestimated=sparse(waitingestimated);

% relocation minutes
Sim.relodist=relodist*P.e;

% trip minutes
Sim.tripdist=tripdist*P.e;

% emissions [ton]
Sim.emissions=(sum(Sim.e/60*P.e,2)')*co2(1:tsim)/10^6;

% pricing info
if isfield(P,'pricing')
    
    distances=Tr(sub2ind(size(Tr),A(:,1),A(:,2)))*P.e; % minutes
    
    Sim.revenues=sum((offeredprices-m.gamma_r).*distances.*chosenmode.*(1-dropped));
    Sim.relocationcosts=sum(relodist)*P.e*m.gamma_r;
    Sim.offeredprices=offeredprices;

    Sim.prices=prices;
    
end



%% create Res struct and save results

% NOTE: need to reorganize results struct

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

if strcmp(P.trlayeralg,'opti') %%P.type<7
    meanqnow=mean(X(n^2+n*P.m*(maxt+2)+1:n^2+n*P.m*(maxt+2)+P.m,  i  ));
else
    meanqnow=mean(q(i,:));
end

fprintf('sim #%d successfully completed - avg soc: %0.2f - total time: %d:%0.2d - \n',dispiter,meanqnow,floor(elapsed/60),round(rem(elapsed,60)));

end






function [fo,fd,Trips]=generateEMD(A,Atimes,T,etsim,TripName,tripday)

    n=size(T,1);
    statsname=['data/temp/tripstats-' TripName '-' num2str(tripday) '-N' num2str(n) '.mat'];
    if exist(statsname,'file')
        load(statsname,'fo','fd','dk');
    else
        [~,fo,fd,dk]=tripstats2(A,Atimes,T);
        save(statsname,'Atimes','fo','fd','dk');
    end


    emdname=['data/temp/emd-' TripName '-' num2str(tripday) '-' num2str(etsim) '.mat'];
    if exist(emdname,'file')
        load(emdname,'dkemd','dkod','dktrip','fk');
    else

        % is a probability distribution of trips available?
        probabilistic=false;

        if probabilistic
            % calculate from known distribution
            error('not implemented');
        else
            % [dkemd,dkod,dktrip,fk]=generatetripdata3(fo,fd,dk,T,etsim);
            [dkemd,dkod,dktrip,fk]=generatetripdataAlt(fo,fd,dk,T,etsim);
        end
        save(emdname,'dkemd','dkod','dktrip','fk');

    end
    
    Trips.dkemd=dkemd;
    Trips.dkod=dkod;
    Trips.dktrip=dktrip;
    Trips.fk=fk;
end

