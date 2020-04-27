% [Res]=GENERALC(P[,extsave,dispiter])
% Run SAEV simulation and relocation/charge otpimization with coordinates.
% Vehicles start at beginning of time step, arrive at end of time step
% 
% CHECK THAT TOTAL DISTANCES TRAVELED BY EV IS THE SAME
% PROVARE TUTTE LE COMBINAZIONI
% CONTROLLA PREVISIONE TRIP: COME � CALCOLATA?
% CHECK P.T: SHOULD BE MINUTES, NOT TIME STEPS
% PUT ALL INFO OF LAYERS INTO STRUCT WITH ALL RELEVANT INFO
% AT WHICH STATIONS VEHICLES ARE CHARGED? (where to put charging stations) 
% FAI DOCUMENTATION
% 
% add model version to Hash?
% 
% what do I still want in Results?
% simulation variables:
% - moving vehicles (p) -> change into more meaningful variable
% decision variables:
% - relocation decisions
% - pickup decisions (associated with each request)
% - charging decision at energy layer
% internals:
% - variables dependent on algorithm
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


%% load trips

% can add secondary trip file (real vs expected/forecasted)
[A,Atimes,ASortInd,AbuckC, ...
    ODistToNode,ONodeID,DDistToNode,DNodeID]=...
    generateGPStrips(P);


%% parameters of simulation

n=size(P.T,1);           % number of nodes
tsim=1440/P.e;          % number of time steps in transport layer
mtsim=tsim/P.beta;      % number of time steps in energy layer
Tr=round(P.T/P.e);      % distance matrix in transport layer steps
ac=round(P.chargekw/P.battery/60*P.e,3);    % charge rate per time step (normalized)
ad=P.consumption/P.battery*P.e;             % discharge rate per time step (normalized)
elep=repelem(P.melep,P.beta);
% maxt=max(max(Tr)); % max distance in time steps

% main variables
q=zeros(tsim,P.m,'double');            % SOC
q(1,:)=P.initialsoc.*ones(1,P.m);      % initial state of charge
e=zeros(tsim,P.m,'double');            % charging
u=zeros(tsim,P.m,'double');            % vehicles in charging stations
u(1,:)=randi(n,1,P.m);               % initial position of vehicles
a=zeros(tsim,n^2,'double');            % number of passengers waiting in stations
p=zeros(P.m,2); % last position
v=zeros(tsim+100,P.m,'double');        % auxiliary variable to assign vehicles to future for trips with passengers
w=zeros(tsim+100,P.m,'double');        % auxiliary variable to assign vehicles to future for relocation
b=zeros(mtsim,n,'double');           % imbalance

queue=zeros(100,1);          % temporary variable to store queued arrivals
% poolid=0;                    % index of pool
waiting=zeros(length(A),1);  % minutes waited for each request
dropped=zeros(length(A),1);  % request is dropped?
traveled=zeros(length(A),1); % trip length (minutes)
pooling=zeros(length(A),1);  % pool ID of each user (if ride shared)


%% trip processing

% generate number of arrivals at each station
statsname=['data/temp/tripstats-' P.tripfile '-' num2str(P.scenarioid) '-N' num2str(n) '.mat'];
if exist(statsname,'file')
    load(statsname,'Atimes','fo','fd','dk');
else
    [Atimes,fo,fd,dk]=tripstats2(A,Atimes,P.T);
    save(statsname,'Atimes','fo','fd','dk');
end

% c=double(convertAtimes(A,Atimes,n,tsim));
% C=reshape(c(1:n^2*tsim),n^2,tsim); % arrivals actual


%% setup energy layer

% generate EMD in case of aggregate energy layer
emdname=['data/temp/emd-' P.tripfile '-' num2str(P.scenarioid) '.mat'];
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

% static variables
Q.minfinalsoc=1; % this only works for optimization horizon of ~24h
Q.storagemax=P.battery*P.m*P.maxsoc; % kWh
Q.maxchargeminute=ac*P.battery;%P.chargekw/60*P.e;  % kWh per time step per vehicle
Q.T=P.mthor;
Q.v2g=P.v2g;
Q.socboost=0;
Q.eta=1;
gridexchange=P.chargekw*P.m*P.beta;
Q.gridimportconstraint=gridexchange;
Q.gridexportconstraint=gridexchange;
Q.selling=1;
Q.cyclingcost=P.cyclingcost;
% P.consumption*P.speedkmh/P.battery/60

% all in units of time steps
Trips.dkmed=dkemd;
Trips.dkod=dkod;
Trips.dktrip=dktrip;
Trips.fk=fk;


%% setup transport layer

zmacro=zeros(4,mtsim+P.mthor); % matrix of optimal control variables for energy layer
relodist=zeros(ceil(tsim/P.tx),1); % distances of relocation

stdepselector=repmat(eye(n),n,1);





















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
    
    if ~strcmp(P.enlayeralg,'no') && rem(i,P.beta)==1 % only cases with energy layer
        
        % current time
        lasttimemacro=cputime;
        
        % index of energy layer
        macroindex=(i-1)/P.beta+1;
        
        if strcmp(P.enlayeralg,'opti')
            
            % ....
            
        elseif strcmp(P.enlayeralg,'aggregate')
            

            if strcmp(P.trlayeralg,'simplified') 
                qnow=q(i,:);
            else
                qnow=(X(n^2+n*P.m*(maxt+2)+1:n^2+n*P.m*(maxt+2)+P.m,  i  ));
            end
            extrasoc=0.15; % extra soc for macro simulation to account for aggregate uncertainty
            actualminsoc=min(P.minsoc+extrasoc,mean(qnow)*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of macro simulation
            Q.storagemin=P.battery*P.m*actualminsoc; % kWh

            % dynamic variables
            % time steps spent traveling during this horizon
            dktripnow=dktrip(macroindex:macroindex+P.mthor-1); 
            Q.einit=sum(qnow)*P.battery;            % total initial energy [kWh]
            Q.etrip=dktripnow*P.consumption;        % energy used per step [kWh] 
            Q.dkav=max(0,P.m*P.beta*P.e-dktripnow); % minutes of availability of cars
            Q.electricityprice=P.melep(macroindex:macroindex+P.mthor-1); 

            E=aevopti11(Q);
            
            % make sure that there is no discharging when V2G is not allowed
            if P.v2g==0
                E.discharging(:)=0;
            end
            
        end
        
        % record time
        S.enlayerCPUtime(macroindex)=cputime-lasttimemacro;
        S.lasttime=cputime;
        
    end
    
    
    
    
    
    
    
    %% transport layer
    
    if strcmp(P.trlayeralg,'simplified') % simplified relocation 
        
        %% relocation
        
        % if it's time for a relocation decision
        if mod(i-1,P.tx)==0
            
            % current relocation number
            kt=(i-1)/P.tx+1;
            
            % number of vehicles at each station
            uv=histc(u(i,:),1:n);
            
            % expected imbalance at stations
            b(kt,:)=uv ...
                -a(i,:)*stdepselector ...  number of passengers waiting at each station
                +sum(fd(i:i+P.ts,:)) ...  expected arrivals between now and now+P.ts
                -sum(fo(i:i+P.ts+P.tr,:)) ...     expected requests between now and now+P.ts+t
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
%             poolingok=[queuepoolingok(queue>0);poolingok];
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
                    uid=(u(i,:)==j);

                    % distance to move
                    distancetomove=Tr(j,A(tripsK,2));

                    % destination station
                    destinations=A(tripsK,2);

                    % sort trips by distance (highest first)
                    [distancetomovesorted,sortid]=sort(distancetomove,'descend');
                    
                    % for each trip
                    for ka=1:length(distancetomovesorted)

                        % trip ID
                        tripID=tripsK(sortid(ka));
                        
                        % find candidate vehicles for the task with enough soc
                        candidates=q(i,:).*uid.*(q(i,:) >= distancetomovesorted(ka)*ad);

                        % vehicle with highest soc goes
                        [ur,ui]=max(candidates);

                        if ur>0
                            
                            % accept request and update vehicle position
                            u(i,ui)=0;
                            uid(ui)=0;
                            v(i+distancetomovesorted(ka),ui)=destinations(sortid(ka));
                            
                            % update distance traveled
%                             Trips.movedist(i)=Trips.movedist(i)+distancetomove;

                        else
                            
                            % check wait time
                            if waiting(tripID)<P.maxwait
                                % if waiting time of this request is less than max. wait time
                                ql=ql+1;  % current trip
                                waiting(tripID)=waiting(tripID)+1;
                                queue(ql)=tripID;
%                                 queuepoolingok(ql)=poolingok(k);
                            else
                                % if max waiting exceeded, request dropped
                                if pooling(tripID)>0
                                    tripsdropped=find(pooling==pooling(tripID));
                                else
                                    tripsdropped=tripID;
                                end
%                                 Trips.rejected(i,A(trips(k),1))=Trips.rejected(i,A(trips(k),1))+length(tripsdropped);
                                dropped(tripsdropped)=1;
                            end
                            
                        end
                    end 
                end
            end
        end

        maxc=Q.dkav*Q.maxchargeminute;
        zmacro(:,macroindex:macroindex+P.mthor-1)=[E.charging./maxc , E.discharging./maxc , maxc , Q.etrip]';
        zmacro(isnan(zmacro))=0;

        v2gallowed=q(i,:)>P.v2gminsoc;
        
        chargevector=(ones(1,P.m)*zmacro(1,macroindex)-v2gallowed*zmacro(2,macroindex))*ac;
        
        % update SOC for vehicles charging
        e(i,:)=min(P.maxsoc,max(P.minsoc,q(i,:)+(u(i,:)>0).*chargevector))-q(i,:);
        
        % update SOC 
        q(i+1,:)=min(P.maxsoc,max(P.minsoc,q(i,:)+e(i,:)-(u(i,:)==0).*ad));
        
        % update idle vehicle positions
        u(i+1,:)=u(i,:)+v(i,:)+w(i,:);
        
        % update waiting passengers
        a(i+1,:)=a(i,:);
        
    end
    
    % record time
    S.trlayerCPUtime(i)=cputime-S.lasttime;
    S.lasttime=cputime;
    
end


%% final calculations

if strcmp(P.trlayeralg,'simplified') 
    
    Sim.d=reshape(a,tsim+1,n,n); % waiting passengers at stations [tsim x n^2]
    Sim.u=u; % final destination of vehicles (station) [tsim x m]
    Sim.q=q; % state of charge 
    Sim.e=e/ac*P.chargekw;
    
    Internals.b=b;
    Internals.v=v;
    Internals.w=w;
    Internals.zmacro=zmacro;
    
%     cfinal=reshape(c(1:n^2*tsim),n^2,tsim);
%     dfinal=a';
    
end


% %% calculate waiting times
% 
% % arrivals [origin/destinationPair  timesteps]
% cfinal=[zeros(n^2,1) cfinal zeros(n^2,1)];
% 
% % waiting at stations [origin/destinationPair  timesteps]
% dfinal=[dfinal zeros(n^2,1)];
% 
% % total arrivals
% totalarrivals=sum(sum(cfinal));
% 
% % cumulative arrivals at each origin/destinationPair
% totalenodi=[0 ; cumsum(sum(cfinal,2))];
% 
% % total waiting minutes
% waitingtimesteps=zeros(totalarrivals,3);
% 
% % for each origin/destinationPair
% for k=1:n^2
%     
%     % find moments with arrivals at origin/destinationPair
%     arrivi=find(cfinal(k,:));
%     
%     % how many arrivals for each moment
%     numeroarrivi=cfinal(k,(cfinal(k,:)>0));
%     
%     % initialize
%     num=1;
%     
%     % for each moment when there are arrivals
%     for k2=1:length(arrivi)
% 
%         for j=1:numeroarrivi(k2)
% 
%             davantiinfila=dfinal(k,:)-(cumsum(cfinal(k,:))-num+1);
%             waitingtimesteps(totalenodi(k)+num,1)=arrivi(k2); % time step of request
%             waitingtimesteps(totalenodi(k)+num,2)=k;         % origin/destination pair
%             waitingtimesteps(totalenodi(k)+num,3)=find((davantiinfila(arrivi(k2):end)<0),1)-1; % waiting time
% 
%             num=num+1;
%             
%         end
%     end
% end
% 
% % waiting time for each passenger in minutes
% waiting=[waitingtimesteps(:,1:2) , waitingtimesteps(:,3)*P.e];
% 
% % moving average [what???]
% binsize=round(10/P.e);
% halfbin=floor(binsize/2);
% waitingprof4=zeros(tsim,1);
% for k2=1:tsim
%     waitingprof4(k2)=mean(    waiting(   logical((waiting(:,1)>=k2-halfbin).*(waiting(:,1)<=k2+halfbin))   ,3)   ,'omitnan');
% end
% waitingprof4(isnan(waitingprof4))=0; % for bins without arrivals, assume 0 waiting



% waiting times
Sim.waiting=waiting;
% Sim.waitingMAV10min=waitingprof4;
% Sim.waitsummary=[  prctile(waiting(:,3),[0,2.5,25,50,75,97.5,100])' ; mean(waiting(:,3)) ; std(waiting(:,3))  ];
% Sim.waitsummaryLgnd=["min";"2.5 pctile";"25 pctile";"50 pctile";"75 pctile";"97.5 pctile";"max";"mean";"st.dev."];

%relocation minutes
Sim.relodist=relodist*P.e;


%% create Res struct and save results

% A(Asortindex,:)=A;
% Atimes(Asortindex,:)=Atimes;



% total cpu time
elapsed=cputime-S.starttime;

% parameters of simulation
Params.Tr=Tr;
Params.elep=elep;
Params.tsim=tsim;

% Res struct generation
Res.Params=Params;
Res.Trips=Trips;
Res.Sim=Sim;
Res.Internals=Internals;
Res.Stats=S;
Res.cputime=elapsed;
Res.cost=(sum(Sim.e/60*P.e,2)')*Params.elep(1:tsim);
Res.peakwait=max(waiting(:,3));
Res.avgwait=mean(waiting(:,3));

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








