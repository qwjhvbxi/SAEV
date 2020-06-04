% %% debug simulation
% 
% P1=cpar('NYC2016-small2','simplified');
% P1.Operations.maxwait=Inf;
% Res1=generalC(P1,2,2)
% P=cpar('NYC2016-small2','opti');
% Res2=generalC(P,2,2)
% 
% return

%% comparison between optimal and simplified

addpath plots

P1=cpar('NYC2016-small','opti');
Res1=generalC(P1,1,2)

P2=cpar('NYC2016-small','simplified');
P2.Operations.maxwait=Inf;
Res2=generalC(P2,1,2)

[P,R]=generateplotline3('NYC2016',[],'Operations.maxwait',[10 Inf]);


%% sensitivity analysis

[P,R]=generateplotline3('NYC2016',[],'Operations.maxwait',[10 Inf]);




%% simulations over 1 year

% use carbon emissions from NYISO 2018, Germany 2019 (more renewables)
% carbon price... 

% P.uinit=Res.Sim.u(end,:)
% P.initialsoc=...;
    
P=cpar('NYC2016');

SOC=zeros(366,P.m);
Uinit=zeros(366,P.m);

SOC(1,:)=ones(1,P.m)*0.7;
uinitfile=['data/uinit/uinit_' P.tripfile '_' num2str(N) '_' num2str(P.K) '.mat'];
if exist(uinitfile,'file')
    load(uinitfile,'uinit');
else
    uinit=findInitialPos(P,A,Atimes,P.K,15);%50
    save(uinitfile,'uinit');
end
Uinit(1,:)=uinit;

% P.Operations.maxwait=Inf;
for k=1:365
    P.tripday=k;
    P.gridday=k;
    P.Operations.initialsoc=SOC(k,:);
    P.Operations.uinit=Uinit(k,:);
    R(k)=generalC(P,1,-k);
    SOC(k+1,:)=R(k).Sim.q(end,:);
    Uinit(k+1,:)=max(1,min(189,R(k).Sim.u(end,:)+sum(R(k).Internals.v(721:end,:))+sum(R(k).Internals.w(721:end,:))));
end


return





%% plot for paper

params
P.trlayeralg='simplified';
P.tx=1;
P.ts=2;
P.tr=5;
P.bmin=0;
P.m=10000;
P.scenario='NY';
load(['data/scenarios/' P.scenario '.mat'],'T');
P.T=round(T);
Res4=general11(P,-1,0);


%%

tsim=720;
x=linspace(0,24,tsim);
x1=linspace(0,24,tsim+1);
x2=linspace(0,24,1440);
xt=0:4:24;

% power exchanged
figure('Units','centimeters','Position',[10,7,10,7])
plot(x,sum(Res4.Sim.e,2)/1000)
ylabel('power (MW)')
hold on
yyaxis right
plot(x,Res4.Params.elep(1:tsim),'r-')
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('electricity price (yen/kWh)')
legend({'power','price'},'Orientation','horizontal')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print -depsc2 figures/power.eps

% waiting 
figure('Units','centimeters','Position',[10,7,10,7])
plot(x,Res4.Sim.waitingMAV10min)
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('average wait time (min.)')
% legend({'power','price'},'Orientation','horizontal','Location','best')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print -depsc2 figures/waiting.eps

% soc 
figure('Units','centimeters','Position',[10,7,10,7])
plot(x1,mean(Res4.Sim.q,2))
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('average SOC')
% legend({'power','price'},'Orientation','horizontal','Location','best')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times');
print -depsc2 figures/soc.eps


%% load trips only

[A,Atimes,ASortInd,AbuckC, ...
    ODistToNode,ONodeID,DDistToNode,DNodeID]=...
    generateGPStrips(P);



