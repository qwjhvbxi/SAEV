% [Res]=GENERALC(P[,extsave,dispiter])
% Run SAEV simulation and relocation/charge otpimization with coordinates.
% Vehicles start at beginning of time step, arrive at end of time step
% 
% CHECK THAT TOTAL DISTANCES TRAVELED BY EV IS THE SAME
% PROVARE TUTTE LE COMBINAZIONI
% CONTROLLA PREVISIONE TRIP: COME è CALCOLATA?
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


%% parameters of simulation

tsim=1440/P.e;    % length of simulation (time steps) (one day)
mtsim=tsim/P.beta; % length of energy layer

ac=round(P.chargekw/P.battery/60*P.e,3);      % charge rate
ad=P.consumption/P.battery*P.e;      % discharge rate (consumptionpermin/battery)
% 
% if iscell(P.scenario)
%     scenarioname=P.scenario{1};
%     tpd=P.scenario{2};
%     scenariotpd=['-' num2str(tpd)];
% else
%     scenarioname=P.scenario;
%     tpd=0;
%     scenariotpd='';
% end
% Tr=round(P.T/P.e);
% n=size(Tr,1);

elep=repelem(P.melep,P.beta);
% maxt=max(max(Tr)); % max distance in time steps

% namearrivals=[DataFolder 'trips_saev/' scenarioname scenariotpd '-' num2str(P.arrivalseed) '.mat'];
% load(namearrivals,'c1','c2')



% 
% % sum arrivals (counted by minute) over one time step
% if P.e>1
%     L=floor(size(c1,3)/2);
%     c1=permute(squeeze(sum(reshape(permute(c1(:,:,1:L*P.e),[3,1,2]),[P.e,L,n,n]),1)),[2,3,1]);
%     c2=permute(squeeze(sum(reshape(permute(c2(:,:,1:L*P.e),[3,1,2]),[P.e,L,n,n]),1)),[2,3,1]);
% end
% 
% % create vectors for optimization (cexpected) and simulation (c)
% c=reshape(c1(:,:,2:tsim+P.thor+1),[n^2*(tsim+P.thor),1]);
% if P.mpcpredict==1
%     cexpected=c;
%     c2=c1;
% else
%     cexpected=reshape(c2(:,:,2:tsim+P.thor+1),[n^2*(tsim+P.thor),1]);
% end


%% setup energy layer

% % generate EMD in case of aggregate energy layer
% emdname=['data/temp/emd-' scenarioname scenariotpd '-' num2str(P.arrivalseed) '-' num2str(mtsim) '.mat'];
% if exist(emdname,'file')
%     load(emdname);
% else
%     % is a probability distribution of trips available?
%     probabilistic=false;
% 
%     if probabilistic
%         % calculate from known distribution
%         nomemap=['data/scenarios/' scenarioname '.mat'];
%         load(nomemap,'OD','TPH');
%         ODreal=OD;
%         for j6=1:24
%             ODreal{j6}(:,1)=OD{j6}(:,1)*TPH(j6)*tpd;
%         end
%         [dkemd,dkod,dktrip,fk]=generatetripdata2(ODreal,Tr,mtsim);
%         if sum(dktrip>P.m*P.beta)>0 % check if there are enough vehicles
%             dktrip((dktrip>P.m*P.beta))
%             warning('need to increase vehicles or decrease TPH')
%         end
%     else
%         % calculate from expected arrivals 
%         for j6=1:24
%             [I,J,V]=find(c2(:,:,j6));
%             ODreal{j6}=[V,I,J];
%         end
%         [dkemd,dkod,dktrip,fk]=generatetripdata2(ODreal,Tr,mtsim);
%     end
%     save(emdname,'dkemd','dkod','dktrip','fk');
% end

% % all in units of time steps
% Trips.dkmed=dkemd;
% Trips.dkod=dkod;
% Trips.dktrip=dktrip;
% Trips.fk=fk;


%% setup transport layer

addpath('../CarSharingModel/functions/stackable');
addpath('../CarSharingModel/functions/');

% load([DataFolder 'trips/' P.tripfile]);
[A,Atimes,Asortindex,AbuckC,...
    OriginDistanceToNode,OriginNodeID,DestinationDistanceToNode,DestinationNodeID]=...
    generateGPStrips(P);

%% initializations

% main variables
q=zeros(tsim,P.m,'double');            % SOC
q(1,:)=P.initialsoc.*ones(1,P.m);      % initial state of charge

e=zeros(tsim,P.m,'double');            % charging
u=zeros(tsim,P.m,'double');            % vehicles in charging stations
u(1,:)=randi(n,1,P.m);               % initial position of vehicles
% a=zeros(tsim,n^2,'double');            % number of passengers waiting in stations

p=zeros(P.m,2); % last position

v=zeros(tsim+100,P.m,'double');        % auxiliary variable to assign vehicles to future for trips with passengers
w=zeros(tsim+100,P.m,'double');        % auxiliary variable to assign vehicles to future for relocation
b=zeros(mtsim,n,'double');           % imbalance

Cexp=reshape(cexpected(1:n^2*tsim),n^2,tsim); % arrivals expected
C=reshape(c(1:n^2*tsim),n^2,tsim); % arrivals actual

C1=reshape(Cexp(1:n*n,:),[n,n,tsim]);
fo=[squeeze(sum(C1,2))' ; zeros(100,n)];
fd=[squeeze(sum(C1,1))' ; zeros(100,n)];

if isfield(P,'tx')
    tx=P.tx;
    ts=P.ts;
    tr=P.tr;
    bmin=P.bmin;
else
    tx=5;
    ts=15;
    tr=15;
    bmin=0;
end

zmacro=zeros(4,mtsim+P.mthor); % matrix of optimal control variables for energy layer
relodist=zeros(ceil(tsim/tx),1); % distances of relocation

stdepselector=repmat(eye(n),n,1);




















%% =============== start of simulation =============== %%

%% variables for progress display and display initializations

S.starttime=cputime;
S.lasttime=S.starttime;
S.trlayerCPUtime=zeros(tsim,1);
S.enlayerCPUtime=zeros(mtsim,1);
if P.v2g==1
    v2gflag='+v2g';
else    
    v2gflag='';
end

% set up time variables for cputime calculation
comptime=[cputime;zeros(tsim,1)];


%% start of iterations

for i=1:tsim
    
    
	%% display progress    
    
    displayState(i,tsim,dispiter,comptime,40)
%     if dispiter>=0
%         elapsed=cputime-S.starttime;
%         expectedtime=elapsed/(i-1)*(tsim-i);
%         if strcmp(P.trlayeralg,'opti')
%             meanqnow=mean(X(n^2+n*P.m*(maxt+2)+1:n^2+n*P.m*(maxt+2)+P.m,  i  ));
%             totwaiting=0; % should implement
%         else
%             meanqnow=mean(q(i,:));
%             totwaiting=sum(a(i,:));
%         end
%         if parcomp  % compact for clusters
%             if rem(i,P.beta)==1
%                 fprintf('sim #%d - avg soc: %5.2f - wait: %d - time: %d:%0.2d / ETA: %d:%0.2d - progress: %2.1f%% \n',dispiter,meanqnow,totwaiting,floor(elapsed/60),round(rem(elapsed,60)),floor(expectedtime/60),round(rem(expectedtime,60)),i/tsim*100)
%             end
%         else        % verbose
%             clc
%             fprintf('\n\n %d / %d \n\n',i,tsim)
%             fprintf('%s/%s%s - avg soc: %5.2f - wait: %d \n\n',P.enlayeralg,P.trlayeralg,v2gflag,meanqnow,totwaiting)
%             fprintf('elapsed: %d:%0.2d \n\n',floor(elapsed/60),round(rem(elapsed,60)))
%             fprintf('ETA: %d:%0.2d \n\n',floor(expectedtime/60),round(rem(expectedtime,60)))
%         end
%     end
    
    
    %% energy layer
    
    if ~strcmp(P.enlayeralg,'no') && rem(i,P.beta)==1 % only cases with energy layer
        
        % current time
        lasttimemacro=cputime;
        
        % index of energy layer
        macroindex=(i-1)/P.beta+1;
        
        if strcmp(P.enlayeralg,'opti')
            
            % load energy layer matrices
            load(namesimmacro);
            
            % assign current values to energy layer variables
            pmaincurrent=reshape(   X(n^2+1:n^2+n*P.m*(maxt+1),i)  ,  [n*P.m,maxt+1]  ); % current p in main simulation,  [n*m  x  distance]
            pcurrent=[pmaincurrent    zeros(n*P.m,(mmaxt+1)*P.beta-(maxt+1))    ];%zeros(n*P.m,0.5*P.beta+(mmaxt+1)*P.beta-(maxt+1))];
            pmacro=sum(reshape(pcurrent(:,1:(mmaxt+1)*(P.beta)),[n,P.m,mmaxt+1,P.beta]),4);
            
            xmacro=[X(1:n^2,i) ; reshape(pmacro,[n*P.m*(mmaxt+1),1]) ; X(n^2+n*P.m*(maxt+1)+1:end,i)]; % create x for energy layer (only p unchanged)
            
            % calculate energy layer matrices
            bdist=bdis*xmacro   +bdisc  +bdisC*mc(n^2*(macroindex-1)+1:n^2*(macroindex-1+P.mthor));
            beqt=beq*xmacro    +beqc   ;%+beqC*c(n^2*(i-1)+1:n^2*(i+P.thor));
            
            % objective function
            f1=  fx*P.beta      +P.rho1*fu    - fsoc*P.rho3;
            
            % temp solution
            if P.v2g==0
                fqv2g=0;
                % add constraints on V2G to Aeq
            end
            
            % add electricity price to objective function
            f=(f1        +P.rho2*((fq+fqv2g*0.9).*(repelem(P.melep(macroindex:macroindex+P.mthor-1),ctrno+ctrnoplus,1))-fqv2g*P.cyclingcost)     )/P.mthor/100; % cyclingcost*P.battery ??
            
            % energy layer optimization
            [zresmacro,~,~,~]=intlinprog(f,intcon,Adis,bdist,Aeq,beqt,lb,ub,options);
            
            % record constraints
            zmacro(:,macroindex:macroindex+P.mthor-1)=reshape(round(zresmacro,3),[ctrno+ctrnoplus,P.mthor]);
            
            % reload previous matrices
            load(namesim)
            
        elseif strcmp(P.enlayeralg,'aggregate')
            
            % time steps spent traveling during this horizon
            dktripnow=dktrip(macroindex:macroindex+P.mthor-1); 

            if strcmp(P.trlayeralg,'simplified') 
                qnow=q(i,:);
            else
                qnow=(X(n^2+n*P.m*(maxt+2)+1:n^2+n*P.m*(maxt+2)+P.m,  i  ));
            end
            extrasoc=0.15; % extra soc for macro simulation to account for aggregate uncertainty
            actualminsoc=min(P.minsoc+extrasoc,mean(qnow)*0.99); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of macro simulation
            gridexchange=P.chargekw*P.m*P.beta;

            % static variables
            Q.minfinalsoc=1; % this only works for optimization horizon of ~24h
            Q.storagemin=P.battery*P.m*actualminsoc; % kWh
            Q.storagemax=P.battery*P.m*P.maxsoc; % kWh
            Q.maxchargeminute=ac*P.battery;%P.chargekw/60*P.e;  % kWh per time step per vehicle
            Q.T=P.mthor;
            Q.v2g=P.v2g;
            Q.socboost=0;
            Q.eta=1;
            Q.gridimportconstraint=gridexchange;
            Q.gridexportconstraint=gridexchange;
            Q.selling=1;
            Q.cyclingcost=P.cyclingcost;
            % P.consumption*P.speedkmh/P.battery/60

            % dynamic
            Q.einit=sum(qnow)*P.battery;                % total initial energy [kWh]
            Q.etrip=dktripnow*ad*P.battery;    % energy used per step [kWh] 
            Q.dkav=max(0,P.m*P.beta-dktripnow);   % steps of availability of cars
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
        if mod(i-1,tx)==0
            
            % current relocation number
            kt=(i-1)/tx+1;
            
            % number of vehicles at each station
            uv=histc(u(i,:),1:n);
            
            % expected imbalance at stations
            b(kt,:)=uv ...
                -a(i,:)*stdepselector ...  number of passengers waiting at each station
                +sum(fd(i:i+ts,:)) ...  expected arrivals between now and now+ts
                -sum(fo(i:i+ts+tr,:)) ...     expected requests between now and now+ts+t
                +histc(reshape(w(i:i+ts,:),P.m*(ts+1),1),1:n)';% vehicles relocating here between now and now+ts
            
            % identify feeder and receiver stations
            F=min(uv,(b(kt,:)-bmin).*(b(kt,:)>=bmin)); % feeders
            R=(-b(kt,:)+bmin).*(b(kt,:)<bmin); % receivers
            
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
        % trips=find((P.A(:,5)==i);
        
        % generate trip requests for this time step
        Requests=(AbuckC(i)+1:AbuckC(i+1))';
        
        
        
        a(i,:)=a(i,:)+C(:,i)';%cexpected((i-1)*n^2+1:i*n^2);
        
        % if there are passengers waiting
        if sum(a(i,:))>0
            
            % for each station
            for ko=1:n
                
                % number of passengers waiting at station ko
                wko=(a(i,ko:n:n^2)); 
                
                % vehicles at this station
                uid=(u(i,:)==ko);
                
                % if there are passengers waiting at station ko
                if sum(wko)>0
                    
                    % destination station
                    dstn=find(wko);
                    
                    % trip lengths
                    arri=Tr(ko,dstn);
                    
                    % duplicate multiple trips
                    arri=repelem(arri,wko(dstn));
                    dstn=repelem(dstn,wko(dstn));
                    [arris,dstnid]=sort(arri,'descend');
                    
                    for ka=1:length(arris)
                        
                        % find candidate vehicles for the task with enough soc
                        candidates=q(i,:).*uid.*(q(i,:) >= arris(ka)*ad);
                        
                        % vehicle with highest soc goes
                        [ur,ui]=max(candidates);
                        
                        if ur>0
                            k=ko+(dstn(dstnid(ka))-1)*n;%od pair
                            a(i,k)=a(i,k)-1;
                            u(i,ui)=0;
                            uid(ui)=0;
                            v(i+arris(ka),ui)=dstn(dstnid(ka));
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
    
    cfinal=reshape(c(1:n^2*tsim),n^2,tsim);
    dfinal=a';
    
end


%% calculate waiting times

% arrivals [origin/destinationPair  timesteps]
cfinal=[zeros(n^2,1) cfinal zeros(n^2,1)];

% waiting at stations [origin/destinationPair  timesteps]
dfinal=[dfinal zeros(n^2,1)];

% total arrivals
totalarrivals=sum(sum(cfinal));

% cumulative arrivals at each origin/destinationPair
totalenodi=[0 ; cumsum(sum(cfinal,2))];

% total waiting minutes
waitingtimesteps=zeros(totalarrivals,3);

% for each origin/destinationPair
for k=1:n^2
    
    % find moments with arrivals at origin/destinationPair
    arrivi=find(cfinal(k,:));
    
    % how many arrivals for each moment
    numeroarrivi=cfinal(k,(cfinal(k,:)>0));
    
    % initialize
    num=1;
    
    % for each moment when there are arrivals
    for k2=1:length(arrivi)

        for j=1:numeroarrivi(k2)

            davantiinfila=dfinal(k,:)-(cumsum(cfinal(k,:))-num+1);
            waitingtimesteps(totalenodi(k)+num,1)=arrivi(k2); % time step of request
            waitingtimesteps(totalenodi(k)+num,2)=k;         % origin/destination pair
            waitingtimesteps(totalenodi(k)+num,3)=find((davantiinfila(arrivi(k2):end)<0),1)-1; % waiting time

            num=num+1;
            
        end
    end
end

% waiting time for each passenger in minutes
waiting=[waitingtimesteps(:,1:2) , waitingtimesteps(:,3)*P.e];

% moving average [what???]
binsize=round(10/P.e);
halfbin=floor(binsize/2);
waitingprof4=zeros(tsim,1);
for k2=1:tsim
    waitingprof4(k2)=mean(    waiting(   logical((waiting(:,1)>=k2-halfbin).*(waiting(:,1)<=k2+halfbin))   ,3)   ,'omitnan');
end
waitingprof4(isnan(waitingprof4))=0; % for bins without arrivals, assume 0 waiting

% waiting times
Sim.waiting=waiting;
Sim.waitingMAV10min=waitingprof4;
Sim.waitsummary=[  prctile(waiting(:,3),[0,2.5,25,50,75,97.5,100])' ; mean(waiting(:,3)) ; std(waiting(:,3))  ];
Sim.waitsummaryLgnd=["min";"2.5 pctile";"25 pctile";"50 pctile";"75 pctile";"97.5 pctile";"max";"mean";"st.dev."];

%relocation minutes
Sim.relodist=relodist*P.e;


%% create Res struct and save results

% total cpu time
elapsed=cputime-S.starttime;

% parameters of simulation
Params.Tr=Tr;
Params.elep=elep;
Params.tsim=tsim;

% information about trip requests
Trips.c=reshape(c(1:n^2*tsim,:),[n,n,tsim]);
Trips.cexpected=reshape(cexpected(1:n^2*tsim),[n,n,tsim]);
% Trips.totalkm=sum(sum(sum(Trips.c,3).*Tr*P.e/60*P.speedkmh));

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








