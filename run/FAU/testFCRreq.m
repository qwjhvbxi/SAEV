function [Res]=testFCRreq(P,var1,var2,var3)

if nargin==3
    R=var1;
    Connected=R.Internals.sc;
    Soc=R.Sim.q(1:end-1,:);
    plotta=var2;
else
    Connected=var1;
    Soc=var2;
    plotta=var3;
end
tsim=size(Connected,1);
days=round(tsim*P.e/1440);
chargekw=P.Tech.chargekw;
battery=P.Tech.battery;
Contract=P.FCR.contracted;
if numel(Contract)>1
    Contract=repelem(Contract(:),1440,1);
end


%% real time requirements

NumIdle=(sum(Connected,2));
Power=NumIdle*chargekw;

% positive FCR, 0% min soc
Energy=sum(min(chargekw,Soc.*Connected*battery),2)/1000;
MinutesDis=Energy./Contract*60; % kWh 

% positive FCR, 20% min soc
Energy20pct=sum(max(0,min(chargekw,(Soc.*Connected-0.2)*battery)),2)/1000; % MW
MinutesDis20pct=Energy20pct./Contract*60; 

% negative FCR
Storage=sum(min(chargekw,(1-Soc.*Connected)*battery),2)/1000;
MinutesCh=Storage./Contract*60; 

% max positive
MaxDis=min(Power/1000,Energy20pct/15*60);

% max negative
MaxCh=min(Power/1000,Storage/15*60);


Res.MaxDis=MaxDis;
Res.MaxCh=MaxCh;


if nargin>2 && plotta>0

%% plots

DataFolder=setDataFolder();
Format='-depsc2';
Resolution=[];
% Format='-dpng';
% Resolution='-r300';

x=linspace(0,24*days,tsim);
% xt=0:4:24;
linevect={'x-','s-','v-','^-','d-','>--','<--','o--','x--'};

figure('Units','centimeters','Position',[10,7,10,4])
hold on
plot(x,Power/1000,'LineWidth',1)
if numel(Contract)>1
    plot(x,Contract,'Color','k','LineWidth',1,'LineStyle','--')
else
    line([0,24*days],[Contract,Contract],'Color','k','LineWidth',1,'LineStyle','--')
end
xlim([0,24*days])
% xticks(xt)
xlabel('hours')
ylabel('power connected (MW)')
% legend({'m=3000','m=4000','m=5000'})
grid on
box on
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
if plotta>1
print([DataFolder 'figures/FAU/required_power'],Format,Resolution);
end

figure('Units','centimeters','Position',[10,7,10,4])
hold on
% plot(x,MinutesDis)
plot(x,MinutesDis20pct)
plot(x,MinutesCh)
line([0,24*days],[15,15],'Color','k','LineWidth',1,'LineStyle','--')
xlim([0,24*days])
ylim([0,400])
% xticks(xt)
% text(10,100,'positive','FontName','Times')
% text(12,300,'negative','FontName','Times')
% text(3,30,'required','FontName','Times')
legend({'positive','negative','required'})
xlabel('hours')
ylabel('max. duration (min)')
% legend({'m=3000','m=4000','m=5000'})
grid on
box on
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
if plotta>1
print([DataFolder 'figures/FAU/required_energy'],Format,Resolution);
end

figure('Units','centimeters','Position',[10,7,10,7])
hold on
% plot(x,MinutesDis)
plot(x,MaxDis,'LineWidth',1.5)
plot(x,MaxCh,'LineWidth',1.5)
xlim([0,24*days])
% xticks(xt)
% text(10,100,'positive','FontName','Times')
% text(12,300,'negative','FontName','Times')
xlabel('hours')
ylabel('max. contract (MW)')
legend({'positive','negative'})
grid on
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
if plotta>1
print([DataFolder 'figures/FAU/possible_contract'],Format,Resolution);
end

end