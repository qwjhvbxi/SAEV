function [Summary,R]=multiDaySim(Period,P,gridoffset)

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

% launch or retrieve simulations
for j=1:length(Period)
    k=Period(j);
    P.tripday=k;
    P.gridday=k+gridoffset;
    P.Operations.initialsoc=SOC(j,:);
    P.Operations.uinit=Uinit(j,:);
    R(j)=generalC(P,1,-j);
    SOC(j+1,:)=R(j).Sim.q(end,:);
    Uinit(j+1,:)=max(1, min( n , ...
        double(R(j).Sim.u(end,:)) ) );
end





% create Summary
Summary.cost=[R(:).cost]';
Summary.dropped=[R(:).dropped]';
Summary.peakwait=[R(:).peakwait]';
Summary.avgwait=[R(:).avgwait]';
Summary.emissions=zeros(length(Period),1);
Summary.waiting=[];
% Summary.pricekm=zeros(length(Period),1);

totreq=0;
totdropped=0;
totwait=0;
totminutes=0;
for j=1:length(Period)
    Summary.emissions(j)=R(j).Sim.emissions;
%     Summary.pricekm(j)=R(j).cost/
    totreq=totreq+length(R(j).Sim.waiting);
    totdropped=totdropped+full(sum(R(j).Sim.dropped));
    totwait=totwait+full(sum(R(j).Sim.waiting));
    totminutes=totminutes+sum(R(j).Sim.tripdist);
    Summary.waiting=[Summary.waiting;R(j).Sim.waiting];
    
end

Summary.totdropped=totdropped/totreq;
Summary.totavgwait=totwait/totreq;
Summary.totwaitprctile=full(prctile(Summary.waiting,[50 97.5 99 100]));
Summary.totminutes=totminutes;

end



