
% Q:
%               storagemax: 70000
%               storagemin: 21000
%          maxchargeminute: 0.16667
%     gridimportconstraint: 5000
%                        T: 24
%                      eta: 0.9
%                     dkav: [24×1 double]
%                    etrip: [24×1 double]
%                  surplus: [24×1 double]
%                      v2g: 1 / 0
%         electricityprice: [24×1 double]
%
% Res
% E=[Q.einit ; zeros(Q.T,1)];     % total energy stored in each time step [kWh]
% charging=zeros(Q.T,1);        % charged energy in time step [kWh]
% discharging=zeros(Q.T,1);     % discharged energy in time step [kWh]
% gridimport=zeros(Q.T,1);      % imported energy [kWh]
% gridexport=x1(3*Q.T+1:Q.T*4);%*Q.selling;
% startups=reshape(round(   x1(   Q.T*4 + repelem(((0:numgenerators-1)*Q.T*4),1,Q.T)+repmat((Q.T*2+1:Q.T*3) ,1,numgenerators)   )         ),Q.T,numgenerators);
% shutdowns=NaN;
% generation=reshape(round(   x1(   Q.T*4 + repelem(((0:numgenerators-1)*Q.T*4),1,Q.T)+repmat((Q.T*3+1:Q.T*4) ,1,numgenerators)   )         ),Q.T,numgenerators);
% totalcost=x1'*f;
%
% added: v2gminsoc, socboost


function [Res]=aevopti11(Q)


%% default values
if ~isfield(Q,'surplus')
    Q.surplus=zeros(Q.T,1);
end
if ~isfield(Q,'generators')
    Q.generators=[];
end
if ~isfield(Q,'carbonprice')
    Q.carbonprice=0; % yen/kg
end
if ~isfield(Q,'startupinit')
    Q.startupinit=[];
end
if ~isfield(Q,'emissionsGridProfile')
    Q.emissionsGridProfile=zeros(Q.T,1);
end


%% initializations
Res.E=[Q.einit ; zeros(Q.T,1)];     % total energy stored in each time step [kWh]
Res.charging=zeros(Q.T,1);        % charged energy in time step [kWh]
Res.discharging=zeros(Q.T,1);     % discharged energy in time step [kWh]
Res.gridimport=zeros(Q.T,1);      % imported energy [kWh]

etripcday=cumsum(Q.etrip); % energy consumption of trips (kWh)


Aeq=[];
beq=[];

% Q.generators [min max cost/kWh costOn ]


numgenerators=size(Q.generators,1);

if Q.v2g==0
    
    % no v2g
    Aeq=[Aeq ; zeros(Q.T) , diag(ones(Q.T,1))  , zeros(Q.T,Q.T*(2+4*numgenerators))];
    beq=[beq ; zeros(Q.T,1)];
    
end

%%    [ charge   discharge  import export [startups    shutdowns   operative  generation]*numgenerators ]

% min max energy
A=[     -tril(ones(Q.T)) , tril(ones(Q.T))/Q.eta  , zeros(Q.T,Q.T*(2+4*numgenerators))];      % min SOC constraint (e(t))
A=[A;   tril(ones(Q.T)) , -tril(ones(Q.T))/Q.eta  , zeros(Q.T,Q.T*(2+4*numgenerators))];      % max SOC constraint (e(t))

% main equation
A=[A; diag(ones(Q.T,1)) , -diag(ones(Q.T,1))      , -diag(ones(Q.T,1))  ,   diag(ones(Q.T,1)), repmat([zeros(Q.T,Q.T*3), -diag(ones(Q.T,1))],1,numgenerators)];

% generators objective functions
fmgs=[ zeros(4*Q.T,1)]; % cost of startup
fmge=[ zeros(4*Q.T,1)]; % cost of energy generation
fco2=[ zeros(2*Q.T,1)   ; Q.emissionsGridProfile*Q.carbonprice ; -Q.emissionsGridProfile*Q.carbonprice ]; % cost of CO2
fsoc=[ -ones(Q.T,1) ; ones(Q.T,1) ; zeros(Q.T*2+4*numgenerators,1) ];

for i=1:numgenerators
    
    Aeq=[Aeq; zeros(Q.T,Q.T*(4 + (i-1)*4)) ,  -diag(ones(Q.T,1))   , +diag(ones(Q.T,1))  , diag(ones(Q.T,1))-diag(ones(Q.T-1,1),-1) , zeros(Q.T)  , zeros(Q.T,Q.T*4*(numgenerators-i))]; % generator status
    
    fmgs=[ fmgs ; (1-Q.startupinit(i))*Q.generators(i,4) ; ones(Q.T-1,1)*Q.generators(i,4)  ; zeros(Q.T*3,1)  ]; % cost of startup
    
    fmge=[ fmge ; zeros(Q.T*3,1) ; ones(Q.T,1)*Q.generators(i,3) ]; % cost of electricity generation
    
    fco2=[ fco2 ; zeros(Q.T*3,1) ; ones(Q.T,1)*Q.generators(i,5)*Q.carbonprice ]; % cost of carbon
    
    A=[A; zeros(Q.T,Q.T*(4 + (i-1)*4)) ,        zeros(Q.T,Q.T*2), -diag(ones(Q.T,1))*Q.generators(i,2) , diag(ones(Q.T,1))          , zeros(Q.T,Q.T*4*(numgenerators-i))]; % generator production min
    A=[A; zeros(Q.T,Q.T*(4 + (i-1)*4)) ,        zeros(Q.T,Q.T*2), diag(ones(Q.T,1))*Q.generators(i,1) , -diag(ones(Q.T,1))          , zeros(Q.T,Q.T*4*(numgenerators-i))]; % generator production max
    
end


% minimum final soc
A=[A; -ones(1,Q.T)      , ones(1,Q.T)             , zeros(1,Q.T*(2+4*numgenerators))];

b=[-etripcday-Q.storagemin+Q.einit; ...% min soc
    etripcday+Q.storagemax-Q.einit; ...% max soc
    Q.surplus  ; % main equation
    zeros(Q.T*numgenerators,1); % gen
    zeros(Q.T*numgenerators,1); % gen
    -etripcday(end)+Q.einit-Q.minfinalsoc*Q.storagemax]; % min final soc


beq=[beq;zeros(numgenerators*Q.T,1)];


maxgen=0;
if numgenerators>0
    maxgen=max(Q.generators(:,2));
end

lb=[zeros(Q.T*(4+numgenerators*4),1)];
ub=[Q.dkav*Q.maxchargeminute;Q.dkav*Q.maxchargeminute; ...
    ones(Q.T,1)*Q.gridimportconstraint; ...
    ones(Q.T,1)*Q.gridexportconstraint; ...
    repmat([ones(Q.T*3,1) ; ones(Q.T,1)*maxgen],numgenerators,1)];


% minimization function [ charge discharge  imports exports [startups       shutdowns operative generation] ]
fbat=[ ones(Q.T,1)*Q.cyclingcost    ;   zeros(Q.T*(3+4*numgenerators),1)  ];

fgrid=[zeros(2*Q.T,1) ;  Q.electricityprice ; -Q.electricityprice*(Q.selling*.99+0.001) ; zeros(Q.T*(4*numgenerators),1) ];



intcon=Q.T*4+repelem(((0:numgenerators-1)*Q.T*4),1,Q.T*3)+repmat([(1:Q.T*3)    ],1,numgenerators);

try
    options = optimoptions('intlinprog','RelativeGapTolerance',1e-4,'Display','none','MaxTime',50);
catch
    options = optimoptions('intlinprog','TolGapRel',0.0001,'MaxTime',50,'Display','none'); % matlab 2015
end
    
f=fbat+fmge+fmgs+fgrid+fco2+fsoc*Q.socboost*mean(Q.electricityprice);

x1=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);



try
    % vehicle charging
    Res.charging=x1(1:Q.T);
    
catch
    
    figure
    plot(Q.surplus);
    hold on
    plot(Q.etrip);
    
    figure
    plot(cumsum(Q.surplus)-cumsum(Q.etrip)+Q.einit-Q.storagemin)
    hold on
    plot((1:Q.T)'*Q.gridimportconstraint)
    plot(cumsum(Q.surplus)-cumsum(Q.etrip)+Q.einit-Q.storagemin+(1:Q.T)'*Q.gridimportconstraint)
    sum(Q.surplus)-sum(Q.etrip)+Q.einit-Q.storagemin+Q.T*Q.gridimportconstraint
    
    x1
    
end


% stored energy
Res.E(2:end)=  Res.E(1)  +cumsum(x1(1:Q.T)) -cumsum(x1(Q.T+1:2*Q.T))/Q.eta  -etripcday;
Res.discharging=x1(Q.T+1:2*Q.T);

% grid import/export
Res.gridimport=x1(2*Q.T+1:Q.T*3);
Res.gridexport=x1(3*Q.T+1:Q.T*4);%*Q.selling;

Res.startups=reshape(round(   x1(   Q.T*4 + repelem(((0:numgenerators-1)*Q.T*4),1,Q.T)+repmat((Q.T*2+1:Q.T*3) ,1,numgenerators)   )         ),Q.T,numgenerators);
Res.shutdowns=NaN;
Res.generation=reshape(round(   x1(   Q.T*4 + repelem(((0:numgenerators-1)*Q.T*4),1,Q.T)+repmat((Q.T*3+1:Q.T*4) ,1,numgenerators)   )         ),Q.T,numgenerators);

Res.totalcost=x1'*f;


end