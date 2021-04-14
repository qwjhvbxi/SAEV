%% [Res,ResSim]=INTEGRATEDFRAMEWORK(P[,E])
% from Boyaci et al., 2015, simplified for operations only. Launch from
% main folder!
% 
% Output/optimization variables:
% n     [T J]       number of available vehicles (full)
% z     [T T J J]   number of vehicles rented from station to station from t to t -> [w] 
% m     [T T J J]   number of dropped requests from station to station from t to t -> [w] 
% q     [T J]       number of vehicles rented from station at beginning of t (full)  
% q'    [T J]       number of vehicles left to station at end of t (full) 
% b     [T]         number of vehicles still rented (moving) (full)
% e     [T]         number of vehicles being relocated (moving) (full)
% r     [T J J]     number of vehicles relocated from station to station at t (full) 

function [Res]=integratedframework(P,E)

T=P.Ib;            % time steps
J=size(P.T,1);     % stations
DataFolder=getdatafolder();
starttime=cputime;
% FileName=[DataFolder 'out/benchmark/outnew-' num2str(J) 'st-' num2str(P.tripday) '-' num2str(P.O) '-' num2str(P.Ib) '-' num2str(P.K) '.mat'];

%% initialize
totaldaylength=1440; % minutes

% load trip file
if isfield(P,'A')
    
    A=P.A;
    Atimes=P.Atimes;
%     fd=P.fd;
%     fo=P.fo;

    
else
    
    load([DataFolder 'trips/' P.tripfile],'A','Atimes','fd','fo');
    
    if iscell(A) && P.tripday>0
        
        A=A{P.tripday};
        Atimes=Atimes{P.tripday};
%         fd=fd{P.tripday};
%         fo=fo{P.tripday};
        
    end
end


A1=[double(A(:,1:2)) , round(double(Atimes(:,1))/P.tau) , round(double(Atimes(:,2))/P.tau)]; % arrival/dep times in minutes
A1=max(1,A1);  % trips have to start at least at time step 1
A1(A1(:,4)>totaldaylength,:)=[];  % remove trips for next day



%% calculate time step lengths and thresholds
tripstartsorted=sort(A1(:,3));
totaltrips=length(A1);
bucketssize=totaltrips/T;
bucketspos=round(0:bucketssize:totaltrips);
Res.thresholds=[0;tripstartsorted(bucketspos(2:end))];
Res.thresholds(end)=totaldaylength;
if totaldaylength==T
    Res.thresholds=(0:totaldaylength)';
end
A2=A1;
for i=1:T
    A2(logical((A1(:,3)>Res.thresholds(i)).*(A1(:,3)<=Res.thresholds(i+1))),3)=i;
    A2(logical((A1(:,4)>Res.thresholds(i)).*(A1(:,4)<=Res.thresholds(i+1))),4)=i;
end

Res.timeintervals=Res.thresholds(2:end)-Res.thresholds(1:end-1);


%% create variables OD & mapping
% variable z,m construction 
% only possible values: contrained by distance from station to station and
% by existance of trips at a certain time.. non-zero OD variables
[Res.A,~,IC]=unique(A2,'rows');
Res.OD=histc(IC,1:max(IC));
w=length(Res.OD);
totalvar=T*J*3  +w*2  +T*J^2  +T*2; % total variables


%% departure and arrival times for relocation
arrivals=[];
for j1=1:J % origin
    for j2=1:J % destination
        arrivals=[arrivals;Res.thresholds(1:end-1)+max(0.01,P.T(j1,j2))];
    end
end

Res.R=[repmat((1:T)',J^2,1) , zeros(T*J^2,1)];
for i=1:T
    Res.R(logical((arrivals>Res.thresholds(i)).*(arrivals<=Res.thresholds(i+1))),2)=i;
end

% relocation distance
Res.RT=reshape(P.T,J^2,1); % distance jl of relocation


clear arrivals bucketssize bucketspos tripstartsorted A1 A2




%% objective function

% number of dropped requests
fdropped=[zeros(T*J+w,1) ; ones(w,1) ; zeros(T*J*2+T*2+T*J^2,1)];

% distance with passengers (minutes)
distances=P.T(sub2ind(size(P.T),Res.A(:,1),Res.A(:,2)));
fmove=[zeros(T*J,1) ; distances ; zeros(w+T*J*2+T*2+T*J^2,1)];

% relocation time steps
% frelocationts=[zeros(T*J*3+w*2+T,1) ; ones(T+T*J^2,1)];

% relocation distance
frelocation=[zeros(T*J*3+w*2+T*2,1) ; repelem(Res.RT,T,1)];

if nargin<2
    E.Cdropped=5;
    E.speed=20;
    E.Ckm=0.1;
    E.tariff=0.25;
end
f=(  fdropped*E.Cdropped  +  fmove*(E.speed/60*E.Ckm-E.tariff)  +  frelocation*E.speed/60*E.Ckm  );

%% disequality constraints
Adis=[];
bdis=[];


%% equality constraints
Aeq=[];
beq=[];





%% constraint 6 z_jl+m_jl=OD_jl     
Aeq=[Aeq; [ sparse( w, T*J )    speye(w)    speye(w)    sparse( w, T*J*2+T*2+T*J^2 ) ] ];
beq=[beq; Res.OD]; 



%% constraint 9a  z_jl^tu - q_j^t = 0      
% tutti i trip che partono da stazione j al tempo u vanno su q_j^u
mapp1=[];
for i1=1:J  % per ogni stazione di partenza
    for i2=1:T   % per ogni departure time
        mapp1=[mapp1;  sparse((Res.A(:,1)==i1).*(Res.A(:,3)==i2))'  ];
    end
end
Aeq=[Aeq; [ sparse( T*J,T*J )    mapp1    sparse( T*J, w )    -speye(T*J)    sparse( T*J, T*J+T*2+T*J^2 ) ] ];
beq=[beq; zeros(T*J,1)]; 



%% constraint 9b  z_jl^tu - q'_j^u = 0    
% tutti i trip che arrivano a stazione j al tempo u vanno su qp_j^u
mapp2=[];
for i1=1:J
    for i2=1:T
        mapp2=[mapp2;  sparse((Res.A(:,2)==i1).*(Res.A(:,4)==i2))'  ];
    end
end
Aeq=[Aeq; [ sparse( T*J,T*J )    mapp2    sparse( T*J, w+T*J )    -speye(T*J)    sparse( T*J, T*2+T*J^2 ) ] ];
beq=[beq; zeros(T*J,1)]; 




%% constraint 10  -n_j +q_j +r_jl -r_lj  <= 0 
temp2=[];
for a2=1:J % destinazione
    temp2a=[];
    for a1=1:T
        temp2a=[ temp2a;    [  sparse((Res.R((a2-1)*T*J+1:a2*T*J,2)==a1).*(Res.R((a2-1)*T*J+1:a2*T*J,1)==a1))'  ] ];
    end
    temp2=[temp2;  [    sparse(T,(a2-1)*T*J) , ... % other stations
                        temp2a , ...
                        sparse(T,(J-a2)*T*J) ] ];
end
constr10=[   -speye(T*J)    ... -n
                sparse(T*J,w*2)  ...  
                speye(T*J)   ...  q
                sparse(T*J,T*J+T*2)   ... 
                repmat(speye(T*J),1,J)-temp2  ]; % r
Adis=[Adis; constr10 ];
bdis=[bdis; zeros(T*J,1)];





%% constraint 11 n_j^t - n_j^(t+1) -q +q' -sum r_jl^t +sum r_jl^u = 0  
temp6=[];
for a2=1:J % destinazione
    temp6a=[];
    for a1=1:T 
        temp6a=[ temp6a ; sparse((Res.R((a2-1)*T*J+1:a2*T*J,2)==a1))'  ];
    end
    temp6=[temp6; [ sparse(T,(a2-1)*T*J) , temp6a  ,    sparse(T,(J-a2)*T*J)       ] ];
end
temp4=spdiags([ones(T*J,1),-ones(T*J,1)],[0 1],T*J,T*J);
constr11=[  temp4 ,  ...  % n - n 
            sparse(T*J,w*2)    ...
            -speye(T*J)    ...-q
            speye(T*J)    ...q'
            sparse(T*J,T*2)  , ...  
            -repmat(speye(T*J),1,J)+temp6  ...-sum_l r  + 
            ];
for a2=1:J
    constr11(a2*T,:)=0;
end
Aeq=[Aeq; constr11 ]; 
beq=[beq; zeros(T*J,1)]; 



%% constraint 12a z - b = 0   
temp7=sparse([]);
for i=1:T
    temp7=[temp7;sparse((Res.A(:,3)<i).*(Res.A(:,4)>=i))'];
end
Aeq=[Aeq; [ sparse(T,T*J)   temp7     sparse(T,w+T*J*2)    -speye(T)     sparse(T,T+T*J^2) ] ];
beq=[beq; zeros(T,1)];  




%% constraint 12b sum -e + r = 0     
temp10=sparse([]);
for i=1:T
    temp10=[temp10;sparse((Res.R(:,1)<i).*(Res.R(:,2)>=i))'];
end
constr12b=[ sparse(T,T*J*3+w*2+T)  ...
            -speye(T)  ...
            temp10    ];
Aeq=[Aeq; constr12b];
beq=[beq; zeros(T,1)]; 



%% constraint 13    sum n+b+e=V    
Aeq=[Aeq; [  repmat(speye(T),1,J)  sparse(T,w*2+T*J*2)    repmat(speye(T),1,2)   sparse(T,T*J^2)  ]];
beq=[beq; ones(T,1)*P.K]; 



% %% constraint 14     RT r <= WH H    
% if isfinite(P.O)
%     Res.WH=Res.timeintervals; % total minutes in a shift at time t
%     Res.h=ones(T,1)*P.O;  % number of relocation personnel in time t
%     Adis=[Adis; [ sparse(T,totalvar-T*J^2)   kron(round(Res.RT)',speye(T)) ] ];
%     bdis=[bdis; Res.WH.*Res.h  ]; 
% end



%% bounds
lb=zeros(totalvar,1);
ub=[];



%% integer constraints
intcon=1:totalvar;




%% additional constraint : initial position
if isfield(P,'uinit')
    addconstinit1=sparse(1,T);
    addconstinit1(1)=1;
    addconstinit=kron(speye(J),addconstinit1);
    Aeq=[Aeq; [ addconstinit sparse(J,totalvar-T*J)] ];
    beq=[beq; P.uinit'];
end



%% options
try
    % matlab 2017
    options=optimoptions('intlinprog',...
        'RelativeGapTolerance',0.02,...
        'IntegerTolerance',1e-4,...
        'MaxTime',1000,...
        'LPMaxIterations',5e5,...
        'Display','iter');%,'MaxNodes',10000,'Heuristics','none');%,'IntegerTolerance',1e-3);
catch
    % matlab 2015 on server
    options=optimoptions('intlinprog',...
        'TolGapRel',0.02,...
        'MaxTime',50000,...
        'RootLPMaxIter',5e6,...
        'LPMaxIter',5e6,...
        'Display','iter'); 
end



%% launch optimization
X=intlinprog(f,intcon,Adis,bdis,Aeq,beq,lb,ub,options);
% tic
% X=linprog(f,Adis,bdis,Aeq,beq,lb,ub);%,optionslin);
% toc

%% generate back variables
X1=round(X);
Res.w=w;
Sol.X=X1;
Sol.n=reshape(X1(1:T*J),T,J);
Sol.z=X1(T*J+1:T*J+w);
Sol.m=X1(T*J+w+1:T*J+w*2);
Sol.q=reshape(X1(T*J+w*2+1:T*J*2+w*2),T,J);
Sol.qp=reshape(X1(T*J*2+w*2+1:T*J*3+w*2),T,J);
Sol.b=X1(T*J*3+w*2+1:T*J*3+w*2+T);
Sol.e=X1(T*J*3+w*2+T+1:T*J*3+w*2+T*2);
Sol.r=reshape(X1(T*J*3+w*2+T*2+1:end),T,J,J);
Res.Sol=Sol;
Res.dropped=sum(Res.Sol.m)/sum(Res.OD);
Res.cputime=cputime-starttime;

% %% generate simulation
% if isfinite(P.O)
%     P.uinit=Res.Sol.n(1,:);
%     P.Oinit=randi(J,P.O,1);
%     P.reloctype=6;
%     P.r=Res.Sol.r;
%     P.timeintervals=Res.timeintervals;
% 
%     ResSim=general5(P);
% else
%     ResSim=NaN;
% end
% save(FileName,'P','Res','ResSim');















