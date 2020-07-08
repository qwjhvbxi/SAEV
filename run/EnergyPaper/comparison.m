%% comparison between optimal and simplified

gridday=13;

P1=cpar('NYC2016-small','opti');
P1.gridfile='NY_DA_2016';
P1.gridday=gridday;
P1.Operations.initialsoc=0.5;
Res1=generalC(P1,1,2)

P2=cpar('NYC2016-small','simplified');
P2.gridfile='NY_DA_2016';
P2.gridday=gridday;
P2.Operations.maxwait=Inf;
P2.Operations.initialsoc=0.5;
% P2.Tech.cyclingcost=100/4000;
Res2=generalC(P2,1,2)

% addpath plots
% [P,R]=generateplotline3('NYC2016',[],'Operations.maxwait',[10 Inf]);


z=linspace(0,24,721);

figure('Units','centimeters','Position',[10,7,10,7])
hold on

yyaxis right
plot(z(1:end-1),Res1.Params.elep(1:720))
ylabel('electricity price ($/MWh)')

yyaxis left
plot(z(1:end-1),sum(Res1.Sim.e,2)/1000)
plot(z(1:end-1),sum(Res2.Sim.e,2)/1000,'k-')
xlim([0,24])
xlabel('hour')
ylabel('power (MW)')
lgnd=legend({'exact','aggregated','price'},'Orientation','vertical','Location','NorthWest');
lgnd.BoxFace.ColorType='truecoloralpha';
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
% print([DataFolder 'figures/JSER/charging'],'-dpng','-r300');



load([DataFolder 'trips/' P2.tripfile])

Reqs=histc(Atimes(:,1),1:1440);
Waits=accumarray([Atimes(:,1);1440],[Res2.Sim.waiting;0]);
plot(Waits)
plot(Reqs)

Waits10min=sum(reshape(Waits,10,144));
Reqs10min=sum(reshape(Reqs,10,144));


x=linspace(0,24,144);

figure('Units','centimeters','Position',[10,7,10,7])
hold on
line([0,24],[0,0])
plot(x,Waits10min./Reqs10min,'k:')
xlim([0,24])
xlabel('hour')
ylabel('average wait time (min.)')
lgnd=legend({'exact','aggregated'},'Orientation','vertical','Location','NorthWest');
lgnd.BoxFace.ColorType='truecoloralpha';
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
% print([DataFolder 'figures/JSER/charging'],'-dpng','-r300');


%%

Res=Res1;

figure('Units','centimeters','Position',[10,7,10,7])
hold on
yyaxis right
plot(z(1:end-1),Res.Params.elep(1:720))
ylabel('electricity price (yen/kWh)')
yyaxis left
stairs(linspace(0,23.5,48),((Res.Internals.zmacro(1,1:48)-Res.Internals.zmacro(2,1:48)).*Res.Internals.zmacro(3,1:48))'*60/30/1000)
plot(z(1:end-1),sum(Res.Sim.e,2)/1000,'k-')
xlim([0,24])
xlabel('hour')
ylabel('power (MW)')
lgnd=legend({'EL power','TL power','price'},'Orientation','vertical','Location','NorthWest');
lgnd.BoxFace.ColorType='truecoloralpha';
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
% print([DataFolder 'figures/JSER/charging'],'-dpng','-r300');

figure('Units','centimeters','Position',[10,7,10,7])
plot(z(1:end-1),sum(R0.Sim.e,2))
xlabel('daily cost (yen)')
ylabel('number of days')
legend({'unscheduled','scheduled'},'Orientation','vertical','Location','NorthWest')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
% print([DataFolder 'figures/JSER/charging'],'-dpng','-r300');



% %% debug simulation
% 
% P1=cpar('NYC2016-small2','simplified');
% P1.Operations.maxwait=Inf;
% Res1=generalC(P1,2,2)
% P=cpar('NYC2016-small2','opti');
% Res2=generalC(P,2,2)
% 
% return