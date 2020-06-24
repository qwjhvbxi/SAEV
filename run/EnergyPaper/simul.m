%% simulations over long period
% sims in two scenarios (NY and Germany)
% 1. without optimized charging
% 2. with optimized charging
% 3. with optimized charging and carbon pricing
% (6x simulations)

% use carbon emissions from NYISO 2018, Germany 2019 (more renewables),
% carbon price... ?

Period=1:4;

% initializations
P=cpar('NYC2016');
P.Operations.maxwait=20;
P.tripfolder='NYC2016'; % trips in 2019
P.gridfile='Germany_DA_2019';
gridoffset=3; % 2016/01/01: Fri; 2019/01/01: Tue; 
P.m=13000;
P.EnergyLayer.mminsoc=0.35;

% no optimization 
P.enlayeralg='no';
[S1,R1]=multiDaySim(Period,P,gridoffset);

P.enlayeralg='aggregate';
[S2,R2]=multiDaySim(Period,P,gridoffset);

P.carbonprice=50; % [$/ton]
[S3,R3]=multiDaySim(Period,P,gridoffset);


figure
hold on
plot(S1.cost+S1.emissions*50,'x-')
plot(S2.cost+S2.emissions*50,'s-')
plot(S3.cost,'o-')

sum(S2.cost+S2.emissions*50)
sum(S3.cost)

sum(S2.emissions)
sum(S3.emissions)

figure
hold on
plot(S1.avgwait)
plot(S2.avgwait)
plot(S3.avgwait)

figure
hold on
plot(S1.dropped)
plot(S2.dropped)
plot(S3.dropped)

figure
hold on
plot(S1.emissions)
plot(S2.emissions)
plot(S3.emissions)



figure
hold on
plot(sum(R2(2).Sim.e,2))
plot(sum(R3(2).Sim.e,2))
yyaxis right
plot(R3(2).Params.co2)



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



