%% [Res]=aevopti11(Q)
% virtual power plant high level optimization.
% Q has inputs: 
%   minfinalsoc
%   storagemin
%   storagemax
%   maxchargeminute
%   T
%   v2g                 boolean
%   eta
%   selling
%   cyclingcost [$/kWh]
%   einit
%   etrip
%   dkav
%   electricityprice [$/kWh]
% optional:
%   socboost
%   gridimportconstraint
%   gridexportconstraint
%   generators          format: [ min , max , cost/kWh , costOn ]
%   startupinit
%   surplus
%   carbonprice [$/ton]
%   emissionsGridProfile [g/kWh]

function [Res]=aevopti11(Q)


%% default values
if ~isfield(Q,'surplus')
    Q.surplus=zeros(Q.T,1);
end
if ~isfield(Q,'generators')
    Q.generators=[];
end
if ~isfield(Q,'carbonprice')
    Q.carbonprice=0;
end
if ~isfield(Q,'startupinit')
    Q.startupinit=[];
end
if ~isfield(Q,'emissionsGridProfile')
    Q.emissionsGridProfile=zeros(Q.T,1);
end
if ~isfield(Q,'gridimportconstraint')
    Q.gridimportconstraint=Inf;
end
if ~isfield(Q,'gridexportconstraint')
    Q.gridexportconstraint=Inf;
end
if ~isfield(Q,'socboost')
    Q.socboost=0;
end
if ~isfield(Q,'minfinalsoc')
    Q.minfinalsoc=0;
end


%% initializations

% initialize results variables
Res.E=[Q.einit ; zeros(Q.T,1)];   % total energy stored in each time step [kWh]
Res.charging=zeros(Q.T,1);        % charged energy in time step [kWh]
Res.discharging=zeros(Q.T,1);     % discharged energy in time step [kWh]
Res.gridimport=zeros(Q.T,1);      % imported energy [kWh]

% initialize optimization variables
% variables: [ charge   discharge  import export [startups    shutdowns   operative  generation]*numgenerators ]
etripcday=cumsum(Q.etrip); % energy consumption of trips (kWh)
numgenerators=size(Q.generators,1);
Aeq=[];
beq=[]; 


%% create optimization


% min max energy
A=[     -tril(ones(Q.T)) , tril(ones(Q.T))/Q.eta  , zeros(Q.T,Q.T*(2+4*numgenerators))];      % min SOC constraint (e(t))
A=[A;   tril(ones(Q.T)) , -tril(ones(Q.T))/Q.eta  , zeros(Q.T,Q.T*(2+4*numgenerators))];      % max SOC constraint (e(t))

% main equation (TODO: need to add curtailment variable in case of large surplus)
%        charge [kWh]            discharge                   import                  export                              generators        
Aeq=[Aeq; diag(ones(Q.T,1)) , -diag(ones(Q.T,1))      , -diag(ones(Q.T,1))  ,   diag(ones(Q.T,1)), repmat([zeros(Q.T,Q.T*3), -diag(ones(Q.T,1))],1,numgenerators)];
beq=[beq;Q.surplus];

% generators objective functions
fmgs=zeros(4*Q.T,1); % cost of startup
fmge=zeros(4*Q.T,1); % cost of energy generation

% cost of CO2
% NOTE: does it make sense to have negative co2 price for v2g??
% fco2=[ zeros(2*Q.T,1)   ; Q.emissionsGridProfile*Q.carbonprice/10^6 ; -Q.emissionsGridProfile*Q.carbonprice/10^6 ]; % cost of CO2
fco2=[ zeros(2*Q.T,1)   ; Q.emissionsGridProfile*Q.carbonprice/10^6 ; zeros(Q.T,1) ]; 

% soc boost
fsoc=[ -ones(Q.T,1) ; ones(Q.T,1)/Q.eta ; zeros(Q.T*2+4*numgenerators,1) ];

% constraints for generators & cost function constructors
for i=1:numgenerators
    
    Aeq=[Aeq; zeros(Q.T,Q.T*(4 + (i-1)*4)) ,  -diag(ones(Q.T,1))   , +diag(ones(Q.T,1))  , diag(ones(Q.T,1))-diag(ones(Q.T-1,1),-1) , zeros(Q.T)  , zeros(Q.T,Q.T*4*(numgenerators-i))]; % generator status
    
    fmgs=[ fmgs ; (1-Q.startupinit(i))*Q.generators(i,4) ; ones(Q.T-1,1)*Q.generators(i,4)  ; zeros(Q.T*3,1)  ]; % cost of startup
    
    fmge=[ fmge ; zeros(Q.T*3,1) ; ones(Q.T,1)*Q.generators(i,3) ]; % cost of electricity generation
    
    fco2=[ fco2 ; zeros(Q.T*3,1) ; ones(Q.T,1)*Q.generators(i,5)*Q.carbonprice ]; % cost of carbon
    
    A=[A; zeros(Q.T,Q.T*(4 + (i-1)*4)) ,        zeros(Q.T,Q.T*2), -diag(ones(Q.T,1))*Q.generators(i,2) , diag(ones(Q.T,1))          , zeros(Q.T,Q.T*4*(numgenerators-i))]; % generator production min
    A=[A; zeros(Q.T,Q.T*(4 + (i-1)*4)) ,        zeros(Q.T,Q.T*2), diag(ones(Q.T,1))*Q.generators(i,1) , -diag(ones(Q.T,1))          , zeros(Q.T,Q.T*4*(numgenerators-i))]; % generator production max
    
end

beq=[beq;zeros(numgenerators*Q.T,1)];

% additional constraints for no v2g
if Q.v2g==0
    Aeq=[Aeq ; zeros(Q.T) , diag(ones(Q.T,1))  , zeros(Q.T,Q.T*(2+4*numgenerators))];
    beq=[beq ; zeros(Q.T,1)];
end

% min final soc
A=[A; -ones(1,Q.T)      , ones(1,Q.T)             , zeros(1,Q.T*(2+4*numgenerators))];

% vector b for all simulation
b=[-etripcday+Q.einit-Q.storagemin; ...% min soc
    etripcday+Q.storagemax-Q.einit; ...% max soc
    zeros(Q.T*numgenerators,1); % gen
    zeros(Q.T*numgenerators,1); % gen
    -etripcday(end)+Q.einit-Q.minfinalsoc*Q.storagemax]; % min final soc

% max generation from generators
if numgenerators>0
    maxgen=max(Q.generators(:,2));
else
    maxgen=0;
end

% lower and upper bounds on variables
lb=zeros(Q.T*(4+numgenerators*4),1);
ub=[Q.dkav*Q.maxchargeminute; ... % charge 
    Q.dkav*Q.maxchargeminute; ... % discharge
    ones(Q.T,1)*Q.gridimportconstraint; ... % import
    ones(Q.T,1)*Q.gridexportconstraint; ... % export
    repmat([ones(Q.T*3,1) ; ones(Q.T,1)*maxgen],numgenerators,1)]; % generators


% cost of battery cycling
fbat=[ ones(Q.T,1)*Q.cyclingcost    ;   zeros(Q.T*(3+4*numgenerators),1)  ];

% cost of grid exchange
% fgrid=[zeros(2*Q.T,1) ;  Q.electricityprice ; -Q.electricityprice*(Q.selling*.99+0.001) ; zeros(Q.T*(4*numgenerators),1) ];
fgrid=[zeros(2*Q.T,1) ;  Q.electricityprice ; -Q.electricityprice*Q.selling+0.001 ; zeros(Q.T*(4*numgenerators),1) ];

% general cost function
f=fbat+fmge+fmgs+fgrid+fco2+fsoc*Q.socboost;

% integer constraints
intcon=Q.T*4+repelem(((0:numgenerators-1)*Q.T*4),1,Q.T*3)+repmat([(1:Q.T*3)    ],1,numgenerators);

try
    options = optimoptions('intlinprog','RelativeGapTolerance',1e-4,'Display','none','MaxTime',50);
catch
    options = optimoptions('intlinprog','TolGapRel',0.0001,'MaxTime',50,'Display','none'); % matlab 2015
end

% optimization 
[x1,~,Flag]=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);

if Flag>0

    % vehicle charging / discharging
    Res.charging=x1(1:Q.T);
    Res.discharging=x1(Q.T+1:2*Q.T);

    % stored energy
    Res.E(2:end)=  Res.E(1)  +cumsum(x1(1:Q.T)) -cumsum(x1(Q.T+1:2*Q.T))/Q.eta  -etripcday;

    % grid import/export
    Res.gridimport=x1(2*Q.T+1:Q.T*3);
    Res.gridexport=x1(3*Q.T+1:Q.T*4);%*Q.selling;

    % generators results
    Res.startups=reshape(round(   x1(   Q.T*4 + repelem(((0:numgenerators-1)*Q.T*4),1,Q.T)+repmat((Q.T*2+1:Q.T*3) ,1,numgenerators)   )         ),Q.T,numgenerators);
    Res.shutdowns=NaN;
    Res.generation=reshape(round(   x1(   Q.T*4 + repelem(((0:numgenerators-1)*Q.T*4),1,Q.T)+repmat((Q.T*3+1:Q.T*4) ,1,numgenerators)   )         ),Q.T,numgenerators);

    Res.totalcost=x1'*f;

else
    
    Res=[];
    
end


end