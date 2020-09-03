%% [Summary,R]=multiDaySim(Period,P[,GridOffset,ResultsOut])
% launch simulations over multiple days.
% Period is the consecutive days
% P is the parameter struct
% GridOffset is the offset of the griddays over the Period
% ResultsOut is a logical variable to indicate if the full results should be
%   outputted (warning: large memory use)
% Summary is a struct with a summary of the results
% R is the optional full results struct if ResultsOut==true

function [Summary,R]=multiDaySim(Period,P,GridOffset,ResultsOut)

if prod(Period(2:end)-Period(1:end-1)==1)==0
    warning('Days in ''Period'' must be consecutive');
    return
end

if nargin<3
    GridOffset=0;
end

if nargin<4
    ResultsOut=false;
end

if  ResultsOut==false
    R=[];
end

DataFolder=setDataFolder();

% initial parameters
SOC=zeros(length(Period)+1,P.m);
Uinit=zeros(length(Period)+1,P.m);
load([DataFolder 'scenarios/' P.scenario '.mat'],'T');
n=size(T,1);

% load or create start variables
StartFile=[DataFolder 'temp/StartFile_' num2str(n) '-' num2str(P.m) '.mat'];
if exist(StartFile,'file')
    load(StartFile,'StartSoc','StartPos');
else
    StartSoc=ones(1,P.m)*0.7;
    StartPos=randi(n,1,P.m);
    save(StartFile,'StartSoc','StartPos');
end
SOC(1,:)=StartSoc;
Uinit(1,:)=StartPos;

Summary.cost=zeros(length(Period),1);
Summary.dropped=zeros(length(Period),1);
Summary.peakwait=zeros(length(Period),1);
Summary.avgwait=zeros(length(Period),1);
Summary.emissions=zeros(length(Period),1);
Summary.waiting=[];
Summary.soc=[];

totreq=0;
totdropped=0;
totwait=0;
totminutes=0;

% launch or retrieve simulations
for j=1:length(Period)
    k=Period(j);
    P.tripday=k;
    P.gridday=k+GridOffset;
    P.Operations.initialsoc=SOC(j,:);
    P.Operations.uinit=Uinit(j,:);
    Res=generalC(P,1,-j);
    SOC(j+1,:)=Res.Sim.q(end,:);
    Uinit(j+1,:)=max(1, min( n , ...
        double(Res.Sim.u(end,:)) ) );
    
    % create Summary
    Summary.cost(j)=Res.cost;
    Summary.dropped(j)=Res.dropped;
    Summary.peakwait(j)=Res.peakwait;
    Summary.avgwait(j)=Res.avgwait;
    Summary.emissions(j)=Res.Sim.emissions;
    
    totreq=totreq+length(Res.Sim.waiting);
    totdropped=totdropped+full(sum(Res.Sim.dropped));
    totwait=totwait+full(sum(Res.Sim.waiting));
    totminutes=totminutes+sum(Res.Sim.tripdist);
    Summary.waiting=[Summary.waiting;Res.Sim.waiting];
    Summary.soc=[Summary.soc;mean(Res.Sim.q,2)];
    
    if ResultsOut
        R(j)=Res;
%         Pmat(j)=P;
    end
end

Summary.totdropped=totdropped/totreq;
Summary.totavgwait=totwait/totreq;
Summary.totwaitprctile=full(prctile(Summary.waiting,[50 97.5 99 100]));
Summary.totminutes=totminutes;

end



