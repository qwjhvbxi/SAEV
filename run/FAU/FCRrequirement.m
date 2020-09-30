% P=cpar('NYC2016-20clusters');
addpath functions utilities plots
DataFolder=setDataFolder();
P=cpar('Munich_clustered');

P.Operations.maxsoc=0.8;

FCR.limits=[49.8,50.2]; % frequency band
FCR.contracted=5;      % MW
FCR.slowchargeratio=0.8;
FCR.fastchargesoc=0.4;
FCR.aggregatechargeratio=0.5;
P.FCR=FCR;

%% test

P.FCR.filename='frequency_30s2015';
Res1=generalC(P,-1,2);


%% real time requirements

Res=Res1;
P0=P;
Contract=10;
% Format='-depsc2';
% Resolution=[];
Format='-dpng';
Resolution='-r300';

tsim=Res.Params.tsim;
NumIdle=(sum(Res.Internals.sc(1:tsim,:),2));
x=linspace(0,24,tsim);
xt=0:4:24;
Power=NumIdle*P0.Tech.chargekw;
Energy=sum(min(P0.Tech.chargekw,Res.Sim.q(1:Res.Params.tsim,:).*Res.Internals.sc*P0.Tech.battery),2)/1000;
Energy20pct=sum(max(0,min(P0.Tech.chargekw,(Res.Sim.q(1:Res.Params.tsim,:).*Res.Internals.sc-0.2)*P0.Tech.battery)),2)/1000;
Storage=sum(min(P0.Tech.chargekw,(1-Res.Sim.q(1:Res.Params.tsim,:).*Res.Internals.sc)*P0.Tech.battery),2)/1000;
MinutesDis=Energy/Contract*60; % kWh 
MinutesDis20pct=Energy20pct/Contract*60; % kWh 
MinutesCh=Storage/Contract*60; % kWh 


%% plots

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


%% test frequency file

if 0
    testfilename='frequency_test';
    f=repelem(repmat([eye(96)*(FCR.limits(1)-50) eye(96)*(FCR.limits(2)-50)],1,2),15,1)+50;
    save([DataFolder 'grid/' testfilename],'f'); % resolution must be at least 1 minute
end

%% launch test

addpath run/FAU

P.FCR=FCR;
P.FCR.filename='frequency_test';
P.m=3000;
P.gridday=17*4+96; % positive from 18:00
% P.gridday=17*4; % negative from 18:00
% P.gridday=16*4+96; % positive from 17:00
Res1=generalC(P,-1,2);

figure
plot(sum(Res1.Sim.ef,2))

figure
hold on
plot(sum(Res1.Sim.e+double(Res1.Sim.ef),2),'k-')
plot(sum(Res1.Sim.e,2),'--')
plot(sum(Res1.Sim.ef,2),'r:')


[FailMinutes,DeltaPower]=testFCR(P,Res1)

