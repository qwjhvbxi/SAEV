function [Res]=testFCRreq(P,R,plotta)

%% real time requirements

Contract=P.FCR.contracted;
tsim=R.Params.tsim;

NumIdle=(sum(R.Internals.sc(1:tsim,:),2));
Power=NumIdle*P.Tech.chargekw;

% positive FCR, 0% min soc
Energy=sum(min(P.Tech.chargekw,R.Sim.q(1:R.Params.tsim,:).*R.Internals.sc*P.Tech.battery),2)/1000;
MinutesDis=Energy/Contract*60; % kWh 

% positive FCR, 20% min soc
Energy20pct=sum(max(0,min(P.Tech.chargekw,(R.Sim.q(1:R.Params.tsim,:).*R.Internals.sc-0.2)*P.Tech.battery)),2)/1000;
MinutesDis20pct=Energy20pct/Contract*60; % kWh 

% negative FCR
Storage=sum(min(P.Tech.chargekw,(1-R.Sim.q(1:R.Params.tsim,:).*R.Internals.sc)*P.Tech.battery),2)/1000;
MinutesCh=Storage/Contract*60; % kWh 

if nargin>2 && plotta

%% plots

% Format='-depsc2';
% Resolution=[];
Format='-dpng';
Resolution='-r300';

x=linspace(0,24,tsim);
xt=0:4:24;
linevect={'x-','s-','v-','^-','d-','>--','<--','o--','x--'};

figure('Units','centimeters','Position',[10,7,10,7])
plot(x,Power/1000)
line([0,24],[Contract,Contract],'Color','red')
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('power connected (MW)')
% legend({'m=3000','m=4000','m=5000'})
grid on
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
% print([DataFolder 'figures/FAU/required_power'],Format,Resolution);

figure('Units','centimeters','Position',[10,7,10,7])
hold on
% plot(x,MinutesDis)
plot(x,MinutesDis20pct)
plot(x,MinutesCh)
line([0,24],[15,15],'Color','red')
xlim([0,24])
xticks(xt)
text(10,100,'positive','FontName','Times')
text(12,300,'negative','FontName','Times')
text(3,30,'required','FontName','Times')
xlabel('hours')
ylabel('max. duration (min)')
% legend({'m=3000','m=4000','m=5000'})
grid on
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
% print([DataFolder 'figures/FAU/required_energy'],Format,Resolution);

end
