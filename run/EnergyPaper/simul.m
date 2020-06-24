%% simulations over long period

% use carbon emissions from NYISO 2018, Germany 2019 (more renewables),
% carbon price... ?

Period=1:30;

% initializations
P=cpar('NYC2016');
P.Operations.maxwait=20;
P.tripfolder='NYC2016';
P.gridfile='Germany_DA_2019';
P.m=13000;
P.EnergyLayer.mminsoc=0.35;
SOC=zeros(length(Period)+1,P.m);
Uinit=zeros(length(Period)+1,P.m);

% initial parameters
load(['data/scenarios/' P.scenario '.mat'],'T');
n=size(T,1);
SOC(1,:)=ones(1,P.m)*0.7;
Uinit(1,:)=randi(n,1,P.m);

%%

for j=1:length(Period)
    k=Period(j);
    P.tripday=k;
    P.gridday=k;
    P.Operations.initialsoc=SOC(k,:);
    P.Operations.uinit=Uinit(k,:);
    R(k)=generalC(P,1,-k);
    SOC(k+1,:)=R(k).Sim.q(end,:);
    Uinit(k+1,:)=max(1, min( n , ...
        double(R(k).Sim.u(end,:)) + full(sum(R(k).Internals.v(721:end,:))) + full(sum(R(k).Internals.w(721:end,:))) ) );
end


%%

Q=[];
for j=1:length(Period)
    k=Period(j);
    Q=[Q;mean(R(k).Sim.q,2)];
    
end

sum([R.cost])

[R.peakwait]

[R.avgwait]

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



