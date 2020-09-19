
%% initializations

addpath functions utilities
DataFolder=setDataFolder();

Period=18:24;
gridfiles={'NY_DA_2016','Germany_DA_2019','Germany_DA_2019'};
gridoffsets=[0,-11,-11+7*20];    % 2016/01/01: Fri; 2019/01/01: Tue;

P=cpar('NYC2016');
P.Operations.maxwait=Inf; %20
P.e=1;
mvec=[5000 6000 8000 10000];

%% create Pmat

P.gridfile=gridfiles{1};
P.carbonprice=0;

P.enlayeralg='no';
P.m=mvec(1);
Pmat{1}=P;
P.m=mvec(2);
Pmat{2}=P;
P.m=mvec(3);
Pmat{3}=P;
P.m=mvec(4);
Pmat{4}=P;

P.enlayeralg='aggregate';
P.m=mvec(1);
Pmat{5}=P;
P.m=mvec(2);
Pmat{6}=P;
P.m=mvec(3);
Pmat{7}=P;
P.m=mvec(4);
Pmat{8}=P;

% P.gridfile=gridfiles{2};
% 
% P.gridfile=gridfiles{3};

GridOffsetsVec=repelem(gridoffsets,1,8);

%% launch multiday sim

% load or create start variables to avoid saving files in parallel computing
load([DataFolder 'scenarios/' P.scenario],'T');
n=size(T,1);
StartFile=[DataFolder 'temp/StartFile_' num2str(n) '-' num2str(P.m) '.mat'];
if exist(StartFile,'file')
    load(StartFile,'StartSoc','StartPos');
else
    StartSoc=ones(1,P.m)*0.7;
    StartPos=randi(n,1,P.m);
    save(StartFile,'StartSoc','StartPos');
end

parfor k=1:length(Pmat)
    multiDaySim(Period,Pmat{k},GridOffsetsVec(k));
end

for k=1:length(Pmat)
    k
    [S1,~]=multiDaySim(Period,Pmat{k},GridOffsetsVec(k));
    S(k)=S1;
end


%% plots

linevect={'x-','s-','v-','^-','d-','>--','<--','o--','x--'};
Format='-depsc2';
Resolution=[];

figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot(mvec,[S(1:4).totavgwait],linevect{1})
plot(mvec,[S(5:8).totavgwait],linevect{2})
ylabel('mean waiting time (min)')
% yyaxis right
% plot(mvec,max([S(1:4).peakwait]),'bx--')
% plot(mvec,max([S(5:8).peakwait]),'rs--')
% ylabel('peak waiting time (min)')
% xlim([0,24])
xticks(mvec)
xlabel('fleet size')
legend({'on-demand','scheduled'},'Orientation','horizontal')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print([DataFolder 'figures/Energy/sens_fleet_wait'],Format,Resolution);


figure('Units','centimeters','Position',[10,7,5,5])
% figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot(mvec,[S(1:4).totavgwait],linevect{1})
plot(mvec,[S(5:8).totavgwait],linevect{2})
ylabel('mean waiting time (min)')
xticks(mvec)
xticklabels([int2str(mvec'/1000) ones(length(mvec),1)*'k'])
xlabel('fleet size')
grid on
legend({'on-demand','scheduled'},'Orientation','vertical')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print([DataFolder 'figures/Energy/sens_fleet_wait-small'],Format,Resolution);

figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot(mvec,mean([S(1:4).cost]),linevect{1})
% plot(([S(1:4).cost]'),'*')
plot(mvec,mean([S(5:8).cost]),linevect{2})
% plot(([S(5:8).cost]'),'*')
ylabel('cost per day ($)')
% xlim([0,24])
xticks(mvec)
xlabel('fleet size')
legend({'on-demand','scheduled'},'Orientation','horizontal','Location','East')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print([DataFolder 'figures/Energy/sens_fleet_cost'],Format,Resolution);

figure('Units','centimeters','Position',[10,7,5,5])
hold on
plot(mvec,mean([S(1:4).cost])/1000,linevect{1})
plot(mvec,mean([S(5:8).cost])/1000,linevect{2})
ylabel('cost per day (thousand $)')
% xlim([0,24])
grid on
xticks(mvec)
xticklabels([int2str(mvec'/1000) ones(length(mvec),1)*'k'])
xlabel('fleet size')
legend({'on-demand','scheduled'},'Orientation','vertical','Location','East')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print([DataFolder 'figures/Energy/sens_fleet_cost-small'],Format,Resolution);

figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot(mvec,mean([S(1:4).cost]),linevect{1})
plot(mvec,mean([S(5:8).cost]),linevect{2})
ylabel('cost per day ($)')
yyaxis right
plot(mvec,[S(1:4).totavgwait],linevect{1})
plot(mvec,[S(5:8).totavgwait],linevect{2})
ylabel('mean waiting time (min)')
% ylabel('peak waiting time (min)')
% xlim([0,24])
xticks(mvec)
xlabel('fleet size')
legend({'on-demand','scheduled'},'Orientation','horizontal')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print([DataFolder 'figures/Energy/sens_fleet_combined'],Format,Resolution);

% k=2;
% % [S0,R]=multiDaySim(18,Pmat{k},GridOffsetsVec(k),1);
% [S1,~]=multiDaySim(Period,Pmat{k},GridOffsetsVec(k),0);