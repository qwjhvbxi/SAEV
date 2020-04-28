% [Res]=GENERALC(P[,extsave,dispiter])
% Run SAEV simulation and relocation/charge otpimization with coordinates.
% Vehicles start at beginning of time step, arrive at end of time step
% 
% TO DO:
% CHECK THAT TOTAL DISTANCES TRAVELED BY EV IS THE SAME
% PROVARE TUTTE LE COMBINAZIONI
% CONTROLLA PREVISIONE TRIP: COME è CALCOLATA?
% PUT ALL INFO OF LAYERS INTO STRUCT WITH ALL RELEVANT INFO
% AT WHICH STATIONS VEHICLES ARE CHARGED? (where to put charging stations) 
% FAI DOCUMENTATION
% add model version to Hash?
% 
% see also PARAMS

function [Res]=generalC(P,extsave,dispiter)


%% initializations

addpath functions
addpath utilities
addpath('../CarSharingModel/functions/stackable');
addpath('../CarSharingModel/functions/');
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


%% load external files: trips and energy

% Note: can add secondary trip file (real vs expected/forecasted)
[A,Atimes,ASortInd,AbuckC,~,~,~,~]=generateGPStrips(P);

load(['data/' P.eleproftype '.mat'],'u','x');
melep=repmat(repelem(x(:,P.eleprofseed),2/P.e,1),2,1);      % macro elep
clear u x;


%% parameters of simulation

% parameters
n=size(P.T,1);              % number of nodes
tsim=1440/P.e;              % number of time steps in transport layer
mtsim=tsim/P.beta;          % number of time steps in energy layer
Tr=max(1,round(P.T/P.e));   % distance matrix in transport layer steps
ac=round(P.chargekw/P.battery/60*P.e,3);    % charge rate per time step (normalized)
ad=P.consumption/P.battery*P.e;             % discharge rate per time step (normalized)
elep=repelem(melep,P.beta);                 % electricity price in each transport layer time step

% main variables
q=zeros(tsim,P.m,'double');            % SOC
e=zeros(tsim,P.m,'double');            % charging
u=zeros(tsim,P.m,'double');            % vehicles in charging stations
v=zeros(tsim+100,P.m,'double');        % auxiliary variable to assign vehicles to future for trips with passengers
w=zeros(tsim+100,P.m,'double');        % auxiliary variable to assign vehicles to future for relocation
b=zeros(mtsim,n,'double');             % imbalance

% working variables
queue=zeros(100,1);          % temporary variable to store queued arrivals
% poolid=0;                    % index of pool

% results variables
waiting=zeros(length(A),1);  % minutes waited for each request
dropped=zeros(length(A),1);  % request is dropped?
pooling=zeros(length(A),1);  % pool ID of each user (if ride shared)
% traveled=zeros(length(A),1); % trip length (minutes)

% initial states
q(1,:)=P.initialsoc.*ones(1,P.m);      % initial state of charge
u(1,:)=randi(n,1,P.m);                 % initial position of vehicles


%% trip processing

% generate number of arrivals at each station
statsname=['data/temp/tripstats-' P.tripfile '-' num2str(P.scenarioid) '-N' num2str(n) '.mat'];
if exist(statsname,'file')
    load(statsname,'fo','fd','dk');
else
    [Atimes,fo,fd,dk]=tripstats2(A,Atimes,P.T);
    save(statsname,'Atimes','fo','fd','dk');
end


%% setup energy layer

% generate EMD in case of aggregate energy layer
emdname=['data/temp/emd-' P.tripfile '-' num2str(P.scenarioid) '-' num2str(mtsim) '.mat'];
if exist(emdname,'file')
    load(emdname,'dkemd','dkod','dktrip','fk');
else
    
    % is a probability distribution of trips available?
    probabilistic=false;

    if probabilistic
        % calculate from known distribution
        error('not implemented');
    else
        % calculate from expected arrivals
        % dkemd, dkod, dktrip are the number of minutes of travel for
        % relocation, serving trips, and total, respectively, for each
        % energy layer time step. fk
        [dkemd,dkod,dktrip,fk]=generatetripdata3(fo,fd,dk,P.T,mtsim);
    end
    save(emdname,'dkemd','dkod','dktrip','fk');
    
end

% energy layer variable: static values
E.minfinalsoc=1;            % this only works for optimization horizon of ~24h
E.socboost=0;
E.v2g=P.v2g;
E.storagemax=P.battery*P.m*P.maxsoc; % kWh
E.maxchargeminute=ac*P.battery;%P.chargekw/60*P.e;  % kWh per time step per vehicle
E.T=P.mthor;                % number of time steps in energy layer
E.eta=1;
E.selling=1;
E.cyclingcost=P.cyclingcost;
% P.consumption*P.speedkmh/P.battery/60

zmacro=zeros(4,mtsim+P.mthor); % matrix of optimal control variables for energy layer
relodist=zeros(ceil(tsim/P.tx),1); % distances of relocation

% all in units of time steps
Trips.dkmed=dkemd;
Trips.dkod=dkod;
Trips.dktrip=dktrip;
Trips.fk=fk;


%% variables for progress display and display initializations

S.starttime=cputime;
S.lasttime=S.starttime;
S.trlayerCPUtime=zeros(tsim,1);
S.enlayerCPUtime=zeros(mtsim,1);

% set up time variables for cputime calculation
comptime=[cputime;zeros(tsim,1)];


%% start of iterations

for i=1:tsim
    
    
	%% display progress    
    
    displayState(i,tsim,dispiter,comptime,40)
    
    
    %% energy layer
    
    if rem(i,P.beta)==1 % only cases with energy layer
        
        % current time
        lasttimemacro=cputime;
        
        % index of energy layer
        macroindex=(i-1)/P.beta+1;
        
        
        switch P.enlayeralg
            
            case 'no'
                
                % charge as much as possible
                zmacro(1,macroindex)=1;
            
            case 'aggregate'
                
                % dynamic variables
                if strcmp(P.trlayeralg,'simplified') 
                    qnow=q(i,:);
                else
                    qnow=(X(n^2+n*P.m*(maxt+2)+1:n^2+n*P.m*(maxt+2)+P.m,  i  ));
                end
                extrasoc=0.15; % extra soc for energy layer to account for aggregate uncertainty
                actualminsoc=min(P.minsoc+extrasoc,mean(qnow)*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of energy layer
                E.storagemin=P.battery*P.m*actualminsoc; % kWh

                dktripnow=dktrip(macroindex:macroindex+P.mthor-1); % time steps spent traveling during this horizon
                E.einit=sum(qnow)*P.battery;            % total initial energy [kWh]
                E.etrip=dktripnow*P.consumption;        % energy used per step [kWh] 
                E.dkav=max(0,P.m*P.beta*P.e-dktripnow); % minutes of availability of cars
                E.electricityprice=melep(macroindex:macroindex+P.mthor-1); 

                ELayerResults=aevopti11(E);

                % make sure that there is no discharging when V2G is not allowed
                if P.v2g==0
                    ELayerResults.discharging(:)=0;
                end

                maxc=E.dkav*E.maxchargeminute;
                zmacro(:,macroindex:macroindex+P.mthor-1)=[ELayerResults.charging./maxc , ELayerResults.discharging./maxc , maxc , E.etrip]';
                zmacro(isnan(zmacro))=0;

            otherwise
            
                error('energy layer must be either ''no'' or ''aggregate'' ');

        end
        
        % record time
        S.enlayerCPUtime(macroindex)=cputime-lasttimemacro;
        S.lasttime=cputime;
        
    end
    
    
    %% transport layer
    
    switch P.trlayeralg
        
        case 'opti'
            
            % need to re-implement
            error('''opti'' transport layer is under implementation');
        
        case 'simplified'       % simplified relocation 
        
            %% charging variables

            v2gallowed=q(i,:)>P.v2gminsoc;
            chargevector=(ones(1,P.m)*zmacro(1,macroindex)-v2gallowed*zmacro(2,macroindex))*ac;


            %% relocation

            % if it's time for a relocation decision
            if mod(i-1,P.tx)==0

                % current relocation number
                kt=(i-1)/P.tx+1;

                % number of vehicles at each station
                uv=histc(u(i,:),1:n);

                % number of waiting passenger at station
                if sum(queue)>0
                    dw=histc(A(queue(queue>0),1)',1:n);
                else 
                    dw=zeros(1,n);
                end

                % expected imbalance at stations
                b(kt,:)=uv ...
                    -dw ...  number of passengers waiting at each station
                    +sum(fd((i-1)*P.e+1:(i+P.ts)*P.e,:)) ...  expected arrivals between now and now+P.ts
                    -sum(fo((i-1)*P.e+1:(i+P.ts+P.tr)*P.e,:)) ...     expected requests between now and now+P.ts+t
                    +histc(reshape(w(i:i+P.ts,:),P.m*(P.ts+1),1),1:n)';% vehicles relocating here between now and now+P.ts

                % identify feeder and receiver stations
                F=min(uv,(b(kt,:)-P.bmin).*(b(kt,:)>=P.bmin)); % feeders
                R=(-b(kt,:)+P.bmin).*(b(kt,:)<P.bmin); % receivers

                % if there are imbalances and available vehicles
                if sum(R)>0 && sum(F)>0

                    % identify optimal relocation flux
                    x=optimalrelocationfluxes(F,R,Tr);

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

                            % save length of relocation
                            relodist(kt)=relodist(kt)+arris(ka);

                        end
                    end
                end
            end
            
        
            %% trip assignment

            % generate trip requests for this time step
            trips=(AbuckC((i-1)*P.e+1)+1:AbuckC(i*P.e+1))';

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

                    % if there are passengers waiting at station k
                    if ~isempty(tripsK)

                        % vehicles at this station
                        uid=find(u(i,:)==j);

                        % soc of vehicles at this station (sorted by high soc)
                        [qj,usortid]=sort(q(i,uid),'descend');

                        % vehicle ID:
                        % uid(usortid)

                        % destination station
                        destinations=A(tripsK,2);

                        % distance to move
                        distancetomove=Tr(j,destinations);

                        % sort trips by distance (highest first)
                        [distancetomovesorted,sortid]=sort(distancetomove,'descend');

                        % for each trip
                        for ka=1:length(distancetomovesorted)

                            % trip ID
                            tripID=tripsK(sortid(ka));

                            % candidates
                            candidates=(qj>=distancetomovesorted(ka)*ad+P.minsoc);

                            % if there are vehicles with enough soc
                            if sum(candidates)>0

                                % vehicle with highest soc goes
                                usortedi=find(candidates,1);

                                % vehicle id
                                ui=uid(usortid(usortedi));

                                % accept request and update vehicle position
                                u(i,ui)=0;
                                qj(usortedi)=0;
                                v(i+distancetomovesorted(ka),ui)=destinations(sortid(ka));

                            else

                                % check if wait time of this request is less than max. wait time (minutes)
                                if waiting(tripID)<P.maxwait

                                    % increase waiting time (minutes) for this trip
                                    ql=ql+1;  % current trip
                                    waiting(tripID)=waiting(tripID)+P.e;

                                    % add this trip to the queue
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
                        end
                    end
                end
            end


            %% simulation variables update

            % update SOC for vehicles charging
            e(i,:)=min(P.maxsoc,max(P.minsoc,q(i,:)+(u(i,:)>0).*chargevector))-q(i,:);

            % update SOC 
            q(i+1,:)=min(P.maxsoc,max(P.minsoc,q(i,:)+e(i,:)-(u(i,:)==0).*ad));

            % update idle vehicle positions
            u(i+1,:)=u(i,:)+v(i,:)+w(i,:);
            
        otherwise
            
            error('transport layer must be either ''opti'' or ''simplified'' ');

    end
    
    % record time
    S.trlayerCPUtime(i)=cputime-S.lasttime;
    S.lasttime=cputime;
    
end


%% final calculations

if strcmp(P.trlayeralg,'simplified') 
    
    Sim.u=u; % final destination of vehicles (station) [tsim x m]
    Sim.q=q; % state of charge 
    Sim.e=e/ac*P.chargekw;
    
    Internals.b=b;
    Internals.v=v;
    Internals.w=w;
    Internals.zmacro=zmacro;
    
end

% waiting times
Sim.waiting=sparse(reorderVectors(waiting,ASortInd));

% dropped requests
Sim.dropped=sparse(reorderVectors(dropped,ASortInd));

%relocation minutes
Sim.relodist=relodist*P.e;


%% create Res struct and save results

% total cpu time
elapsed=cputime-S.starttime;

% parameters of simulation
Params.Tr=uint8(Tr);
Params.elep=elep;
Params.tsim=tsim;

% Res struct generation
Res.Params=Params;
Res.Trips=Trips;
Res.Sim=Sim;
Res.Internals=Internals;
Res.Stats=S;
Res.cputime=elapsed;
Res.cost=(sum(Sim.e/60*P.e,2)')*elep(1:tsim);
Res.dropped=sum(dropped)/length(A);
Res.peakwait=max(waiting);
Res.avgwait=mean(waiting);

% save results
if extsave>0
    save(simname,'Res','P');
end


%% end display

if ~strcmp(P.trlayeralg,'simplified') %%P.type<7
    meanqnow=mean(X(n^2+n*P.m*(maxt+2)+1:n^2+n*P.m*(maxt+2)+P.m,  i  ));
else
    meanqnow=mean(q(i,:));
end

fprintf('sim #%d successfully completed - avg soc: %0.2f - total time: %d:%0.2d - \n',dispiter,meanqnow,floor(elapsed/60),round(rem(elapsed,60)));


end








