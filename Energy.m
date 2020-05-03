%% launch simulation

P=cpar('NYC2018');
P.scenarioid=1;
P.e=2;
P.Operations.maxwait=0;
% P.enlayeralg='no';
% P.consumption=0;
Res2=generalC(P,-1,2)


%% launch over several days

% P.uinit=Res.u(end);
% P.initialsoc=...;
    

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
set(gca,...
...'Units','normalized',...
...'Position',[.15 .2 .75 .7],...
'FontUnits','points',...
'FontWeight','normal',...
...'FontSize',9,...
'FontName','Times')
print -depsc2 figures/power.eps

% waiting 
figure('Units','centimeters','Position',[10,7,10,7])
plot(x,Res4.Sim.waitingMAV10min)
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('average wait time (min.)')
% legend({'power','price'},'Orientation','horizontal','Location','best')
set(gca,...
...% 'Units','normalized',...
...% 'Position',[.15 .2 .75 .7],...
'FontUnits','points',...
'FontWeight','normal',...
...'FontSize',9,...
'FontName','Times')
print -depsc2 figures/waiting.eps

% soc 
figure('Units','centimeters','Position',[10,7,10,7])
plot(x1,mean(Res4.Sim.q,2))
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('average SOC')
% legend({'power','price'},'Orientation','horizontal','Location','best')
set(gca,...
...'Units','normalized',...
...'Position',[.15 .2 .75 .7],...
'FontUnits','points',...
'FontWeight','normal',...
...'FontSize',9,...
'FontName','Times')
print -depsc2 figures/soc.eps


%% load trips only

[A,Atimes,ASortInd,AbuckC, ...
    ODistToNode,ONodeID,DDistToNode,DNodeID]=...
    generateGPStrips(P);



