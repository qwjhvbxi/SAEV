% [Res]=GENERALMILP(P[,extsave,dispiter])
% Run SAEV simulation and relocation/charge oprimization.
% Vehicles start at beginning of time step, arrive at end of time step.
% 
% see also CPAR

function [Res]=generalmilp(P,extsave,dispiter)


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

load([DataFolder 'scenarios/' P.scenario '.mat'],'T');

% NOTE: can add secondary trip file (real vs expected/forecasted)
% [A,Atimes,~,~]=loadTrips(P);
[A,Atimes,AbuckC,~]=loadTrips(P);
% AbuckC=AbuckC(1:P.Sim.e:end);
if ~P.mpcpredict
    Pb.tripfolder=[P.tripfolder '-2'];
    Pb.tripday=P.tripday;
    [A2,Atimes2,~,~]=loadTrips(Pb);
end

% NOTE: should generalize vector length for cases with different beta, e,
% etc. Also: change names of variables
% elep is in $/MWh
load([DataFolder 'eleprices/' P.gridfile '.mat'],'x','y');


%% parameters of simulation

% parameters
n=size(T,1);              % number of nodes
r=AbuckC(1441);           % number of requests
tsim=1440/P.e;            % number of time steps
etsim=floor(1440/P.beta); % number of charging decisions
Tr=max(1,round(T/P.e));   % distance matrix in steps
ac=round(P.Tech.chargekw/P.Tech.battery/60*P.e,3);    % charge rate per time step (normalized)
ad=P.Tech.consumption/P.Tech.battery*P.e;             % discharge rate per time step (normalized)
mthor=round(P.EnergyLayer.mthor/P.beta);

% main variables
q=zeros(tsim,P.m,'double');            % SOC
e=zeros(tsim,P.m,'double');            % charging
u=zeros(tsim,P.m,'double');            % vehicles in charging stations

% results variables
waiting=zeros(r,1);  % minutes waited for each request
dropped=zeros(r,1);  % request is dropped?
chosenmode=false(r,1);% which mode is chosen?
relodist=zeros(ceil(tsim),1); % distances of relocation (at moment of decision)
tripdist=zeros(ceil(tsim),1); % distances of trips (at moment of acceptance)
waitingestimated=zeros(r,1);  % estimated minutes to wait for each request

% initial states
q(1,:)=P.Operations.initialsoc.*ones(1,P.m);      % initial state of charge
if isfield(P.Operations,'uinit')
    u(1,:)=P.Operations.uinit;
else
    u(1,:)=randi(n,1,P.m);                 % initial position of vehicles
end

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


%% setup energy layer

if strcmp(P.enlayeralg,'aggregate') 

    % generate aggregate trip statistics
    EMDFileName=[TripName '-' num2str(P.tripday)];
    [Trips]=computeemd(A,Atimes,T,P.beta,EMDFileName);

%     % append values for next day
%     if isfield(P,'tripfolder')
%         P2=P;
%         P2.tripday=P.tripday+1;
%         [A2,Atimes2,~,~]=loadTrips(P2);
%         EMDFileName=[TripName '-' num2str(P2.tripday)];
%         [Trips2]=computeemd(A2,Atimes2,T,etsim,EMDFileName);
%         Trips.dktrip=[Trips.dktrip(1:48,:) ; Trips2.dktrip(1:48,:)];
%     end
    
    % energy layer variable: static values
    E.v2g=P.Operations.v2g; % use V2G?
    E.eta=1;                % 
    E.selling=1;            % can sell to the grid?
%     E.minfinalsoc=0.9;      % final SOC. This only works for optimization horizon of ~24h
    E.socboost=10000;
    E.T=mthor;% number of time steps in energy layer
    E.cyclingcost=P.Tech.cyclingcost;                       % battery cycling cost [$/kWh]
    E.storagemax=P.Tech.battery*P.m*P.Operations.maxsoc;    % max total energy in batteries [kWh]
    E.maxchargeminute=P.Tech.chargekw/60;                   % energy exchangeable per minute per vehicle [kWh]
    E.carbonprice=P.carbonprice;                            % carbon price [$ per kg]

    % matrix of optimal control variables for energy layer
    zmacro=zeros(4,etsim+mthor); 
    
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
    cname=[DataFolder 'temp/c-' TripName '-' num2str(P.tripday) '.mat'];
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
    
    % create secondary file if needed
    if P.mpcpredict==true
        cexpected=c;
    else
        
        c2=double(convertAtimes(A2,Atimes2,n));
        c21=cat(3,c2,zeros(n,n,tsim)); % add padding
        
        % sum arrivals (counted by minute) over one time step
        if P.e>1
            L=floor(size(c21,3)/2);
            c21=permute(squeeze(sum(reshape(permute(c21(:,:,1:L*P.e),[3,1,2]),[P.e,L,n,n]),1)),[2,3,1]);
            % NOTE: implementation for stochastic arrivals needed:
            % c2=permute(squeeze(sum(reshape(permute(c2(:,:,1:L*P.e),[3,1,2]),[P.e,L,n,n]),1)),[2,3,1]);
        end
        cexpected=reshape(c21(:,:,1:tsim+P.TransportLayer.thor),[n^2*(tsim+P.TransportLayer.thor),1]);

%         % NOTE: need to implement case with imperfect predictions
%         error('mpcpredict==false not implemented'); 
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
    % EMaxHorizon=mthor;
    
end


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
    
    if rem(i,P.beta/P.e)==1 % only cases with energy layer
        
        % current time
        lasttimemacro=cputime;
        
        % index of energy layer
        t=(i-1)/(P.beta/P.e)+1;
                
        switch P.enlayeralg
            
            case 'no'
                
                % charge as much as possible
                zmacro(1,t)=1;
            
            case 'aggregate'
                
                % dynamic variables
                actualminsoc=min(P.Operations.minsoc+P.EnergyLayer.extrasoc,mean(q(i,:))*0.9); % soft minsoc: to avoid violating contraints in cases where current soc is lower than minsoc of energy layer
                E.storagemin=P.Tech.battery*P.m*actualminsoc; % kWh

                dktripnow=Trips.dktrip(t:t+mthor-1); % time steps spent traveling during this horizon
                E.einit=sum(q(i,:))*P.Tech.battery;            % total initial energy [kWh]
                E.etrip=dktripnow*P.Tech.consumption;        % energy used per step [kWh] 
                E.dkav=max(0,P.m*P.beta-dktripnow); % minutes of availability of cars
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
    
%     v2gallowed=q(i,:)>P.Operations.v2gminsoc;
%     chargevector=(ones(1,P.m)*zmacro(1,t)-v2gallowed*zmacro(2,t))*ac;
    
    
    %% transport layer
    
    switch P.trlayeralg
        
        case 'opti'
            
            % adjust values of known parameters (vector b): multiply by current Param.x and add static values
            bdist=bdis*X(1:varno,i)  +  bdisc  + ...
                bdisC * [  c(n^2*(i-1)+1:n^2*i)  ; cexpected(n^2*(i)+1:n^2*(i+P.TransportLayer.thor-1))]; % bdist must be positive
            beqt=beq*X(1:varno,i)    +beqc   ;%+beqC*c(n^2*(i-1)+1:n^2*(i+P.TransportLayer.thor));
                       
            % determine which vehicles can discharge to grid
%             v2gallowed=[ones(P.m,1);(q(i,:)'>P.Operations.v2gminsoc)]*ones(1,EMaxHorizon);
            v2gallowed=[ones(P.m,1);(q(i,:)'>P.Operations.v2gminsoc)];
            
            % create charge vector
%             chargevector=repmat(     reshape(repelem(zmacro(1:2,t:t+EMaxHorizon-1),P.m,1).*v2gallowed,EMaxHorizon*2*P.m,1),P.beta,1);
            chargevector=repmat( reshape(repelem(zmacro(1:2,t),P.m,1).*v2gallowed,2*P.m,1),P.TransportLayer.thor,1);

            % apply charging constraints
%             ub(ubChargingSelector)=chargevector(1:P.TransportLayer.thor*P.m*2)*ac;
            ub(ubChargingSelector)=chargevector*ac;
        
            % transport layer optimization
            zres=intlinprog(f,intcon,Adis,bdist,Aeq,beqt,lb,ub,options);
            
            % optimal control variables in this time step (round binary values)
            z(:,i)=[round(zres(1:P.m*n^2*2));zres(P.m*n^2*2+1:ctrno)];
            
            % calculate new Param.x (using first time step of solution)
            X(1:varno,i+1)=round(Aopti*X(1:varno,i)+Bopti*z(:,i)+B(:,i),4);
            
            % calculate new SOC
            q(i+1,:)=(X(n^2+n*P.m*(maxt+2)+1:n^2+n*P.m*(maxt+2)+P.m,  i+1  ));
            
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
Res.dropped=sum(dropped)/r;
Res.peakwait=max(waiting);
Res.avgwait=mean(waiting);

% save results
if extsave>0
    save(simname,'Res','P');
end


%% end display

meanqnow=mean(q(i,:));

fprintf('sim #%d successfully completed - avg soc: %0.2f - total time: %d:%0.2d - \n',dispiter,meanqnow,floor(elapsed/60),round(rem(elapsed,60)));

end


%% [A,Atimes]=CONVERTATIMES(passengers) 
% from matrix form of OD to list form
% 
% [passengers]=CONVERTATIMES(A,Atimes,n,tsim) 
% from list form to matrix form
% 

function [a,b]=convertAtimes(c,d,e)

if nargin<2
    
    % from list form to matrix form
    
    passengers=c;
    
    A=NaN;
    Atimes=NaN;
    
    a=A;
    b=Atimes;
    
    warning('not implemented');
    return
    
    
else
    
    % from matrix form of OD to list form
    
    A=c;
    Atimes=d;
    n=e; % number of nodes
    tsim=1440; % number of time steps
    b=NaN;
    
    a=zeros(n,n,tsim,'uint16');

%     for i=1:length(A)
%         a(A(i,1),A(i,2),Atimes(i,1))=a(A(i,1),A(i,2),Atimes(i,1))+1;
%     end
    
    for i=1:tsim
        kt=(Atimes(:,1)==i);
        
        A2=sub2ind([n,n],A(kt,1),A(kt,2));
        a(:,:,i)=reshape(accumarray(A2,1,[n^2,1]),n,n);
    end
    
end

end