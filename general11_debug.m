%% debug file

clear

% default
params

P.trlayeralg='opti';
Res1=general11(P,1,0);

P.trlayeralg='simplified';
P.tx=1;
P.ts=2;
P.tr=5;
P.bmin=0;

Res2=general11(P,2,-1);







figure
subplot(2,1,1)
plot(Res1.Sim.q)
subplot(2,1,2)
plot(Res2.Sim.q)

limmin=-(sum(sum(Res1.Sim.e<0))+sum(sum(Res2.Sim.e<0))>0)*P.chargekw*1.1;

figure
subplot(2,1,1)
plot(Res1.Sim.e)
ylim([limmin,P.chargekw*1.1])
subplot(2,1,2)
plot(Res2.Sim.e)
ylim([limmin,P.chargekw*1.1])

figure
plot(Res1.Sim.waitingMAV10min)
hold on
plot(Res2.Sim.waitingMAV10min)

%% new plots

tsim=720;
x=linspace(0,24,tsim);
x1=linspace(0,24,tsim+1);
x2=linspace(0,24,1440);
xt=0:4:24;

% power exchanged
figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot(x,sum(Res1.Sim.e,2)/1000)
plot(x,sum(Res2.Sim.e,2)/1000)
ylabel('power (MW)')
% yyaxis right
% plot(x,Res1.Params.elep(1:tsim),'r-')
xlim([0,24])
xticks(xt)
xlabel('hours')
% ylabel('electricity price (yen/kWh)')
legend({'exact model','aggregated'},'Orientation','horizontal','Location','best')
set(gca,...
...'Units','normalized',...
...'Position',[.15 .2 .75 .7],...
'FontUnits','points',...
'FontWeight','normal',...
...'FontSize',9,...
'FontName','Times')
print -depsc2 figures/C_power.eps

% waiting 
figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot(x,Res1.Sim.waitingMAV10min)
plot(x,Res2.Sim.waitingMAV10min)
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('average wait time (min.)')
legend({'exact model','aggregated'},'Orientation','vertical','Location','West')
set(gca,...
...% 'Units','normalized',...
...% 'Position',[.15 .2 .75 .7],...
'FontUnits','points',...
'FontWeight','normal',...
...'FontSize',9,...
'FontName','Times')
print -depsc2 figures/C_waiting.eps

% soc 
figure('Units','centimeters','Position',[10,7,10,7])
plot(x1,mean(Res1.Sim.q,2))
hold on
plot(x1,mean(Res2.Sim.q,2))
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('average SOC')
legend({'exact model','aggregated'},'Orientation','vertical','Location','SouthWest')
set(gca,...
...'Units','normalized',...
...'Position',[.15 .2 .75 .7],...
'FontUnits','points',...
'FontWeight','normal',...
...'FontSize',9,...
'FontName','Times')
print -depsc2 figures/C_soc.eps


