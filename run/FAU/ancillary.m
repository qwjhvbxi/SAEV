%% launch simulations for MT-ITS 2021 paper 2
% 
% 

%% generate variables

addpath functions utilities plots run/FAU
DataFolder=setDataFolder();
P=cpar('Munich_clustered_week');

P.gridfile='Germany_DA_2015';
P.Operations.maxsoc=0.9;
P.Tech.efficiency=0.9;

% FCR provision parameters
FCR.filename='frequency_30s2015';
FCR.limits=[49.8,50.2]; % frequency band
FCR.contracted=5;      % MW
FCR.slowchargeratio=0.8;
FCR.fastchargesoc=0.4;
FCR.aggregatechargeratio=0.5;

% comparison with set point without FCR 
FCR2=FCR;
FCR2.contracted=0;      % MW
FCR2.slowchargeratio=1;
FCR2.fastchargesoc=1;
FCR2.aggregatechargeratio=1;


%% launching sims over 1 week

% fleet: 5000
P.m=5000;
Pmat{1}=P;
Pmat{1}.enlayeralg='night';

Pmat{2}=P;
Pmat{2}.FCR=FCR2;

Pmat{3}=P;
Pmat{3}.FCR=FCR;

% smaller fleet
P.m=4500;
% P.m=4750;
Pmat{4}=P;
Pmat{4}.enlayeralg='night';

Pmat{5}=P;
Pmat{5}.FCR=FCR2;

Pmat{6}=P;
Pmat{6}.FCR=FCR;
Pmat{6}.FCR.contracted=1;

% % fleet: 6000
% P.m=6000;
% Pmat{7}=P;
% Pmat{7}.enlayeralg='no';
% 
% Pmat{8}=P;
% Pmat{8}.FCR=FCR2;
% 
% Pmat{9}=P;
% Pmat{9}.FCR=FCR;
% Pmat{9}.FCR.contracted=10;


Pmat{7}

%%

Period=1:7;
GridOffset=12-1; % Monday 12 Jan - Sunday 18 Jan
ResultsOut=false;

% StartDay=26+7*5;
% Settimana=(StartDay-1)*24+1:(StartDay+6)*24; % Monday 12 Jan - Sunday 18 Jan
% plot(DA_DE_2015(Settimana));

parfor i=1:length(Pmat)
[~,~]=multiDaySim(Period,Pmat{i},GridOffset,ResultsOut);
end

for i=1:length(Pmat)
[Sres,~]=multiDaySim(Period,Pmat{i},GridOffset,ResultsOut);
S(i)=Sres;
end

% results over 1 week / different seasons
% effect of different fleet sizes, power connection on charging costs, FCR provision
% effect of large FCR request (max for 15 minutes) on passenger transport quality
% synergies of larger fleet size in terms of better transport quality and higher contracted FCR
% loss of charging efficiency with enforcement of constant set point

% Res1=generalC(P1,1,2)
% Res2=generalC(P2,1,2)
% 
% addpath plots
% plotta(Res1,'power')
% plotta(Res2,'power')

%% case with variable FCR contracts

FCR3=FCR;
FCR3.contracted=0;

Period=1:7;
GridOffset=12-1; % Monday 12 Jan - Sunday 18 Jan
ResultsOut=false;

P.m=5000;
Pmat2{1}=P;
Pmat2{1}.FCR=FCR3;

P.m=4500;
Pmat2{2}=P;
Pmat2{2}.FCR=FCR3;


parfor i=1:length(Pmat2)
[~,~]=multiDaySim(Period,Pmat2{i},GridOffset,ResultsOut);
end

for i=1:length(Pmat2)
[Sres,~]=multiDaySim(Period,Pmat2{i},GridOffset,ResultsOut);
S1(i)=Sres;
end

R1=testFCRreq(Pmat2{1},S1(1).sc,S1(1).socv,0);
VariableContract(1,:)=floor(min(min(reshape(R1.MaxDis,1440,7)),min(reshape(R1.MaxCh,1440,7))));

R2=testFCRreq(Pmat2{2},S1(2).sc,S1(2).socv,0);
VariableContract(2,:)=floor(min(min(reshape(R2.MaxDis,1440,7)),min(reshape(R2.MaxCh,1440,7))));

% R1=testFCRreq(Pmat{2},S(2).sc,S(2).socv,0);
% VariableContract(1,:)=floor(min(min(reshape(R1.MaxDis,1440,7)),min(reshape(R1.MaxCh,1440,7))));
% R2=testFCRreq(Pmat{5},S(4).sc,S(4).socv,0);
% VariableContract(2,:)=floor(min(min(reshape(R2.MaxDis,1440,7)),min(reshape(R2.MaxCh,1440,7))));

% prova
Pmat2{2}.EnergyLayer.extrasoc=0.3;

parfor i=1:length(Pmat2)
[~,~]=multiDaySim(Period,Pmat2{i},GridOffset,ResultsOut,'FCR.contracted',VariableContract(i,:));
end

for i=1:length(Pmat2)
[Sres,~]=multiDaySim(Period,Pmat2{i},GridOffset,ResultsOut,'FCR.contracted',VariableContract(i,:));
S2(i)=Sres;
end

Ptest=Pmat2{1};
Ptest.FCR.contracted=VariableContract(1,:);
testFCRreq(Ptest,S2(1).sc,S2(1).socv,1);

Ptest=Pmat2{2};
Ptest.FCR.contracted=VariableContract(2,:);
testFCRreq(Ptest,S2(2).sc,S2(2).socv,1);






%% frequency plot

DataFolder=setDataFolder();
load([DataFolder 'grid/frequency_30s2015.mat'],'f');
Format='-depsc2';
Resolution=[];

figure('Units','centimeters','Position',[10,7,10,7])
v=linspace(0,24*7,7*size(f,1));
plot(v,reshape(f(:,12:18),7*size(f,1),1))
xlim([0,24*7])
ylim([49.9,50.1])
xlabel('hour');
ylabel('frequency (Hz)')
box on
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print([DataFolder 'figures/FAU/frequency'],Format,Resolution);

%% electricity prices

load([DataFolder 'eleprices/Germany_DA_2015.mat'],'x');
Format='-depsc2';
Resolution=[];

figure('Units','centimeters','Position',[10,7,10,7])
v=linspace(0,24*7,7*size(x,1));
plot(v,reshape(x(:,12:18),7*size(x,1),1))
xlim([0,24*7])
xlabel('hour');
ylabel('prices (EUR/MWh)')
box on
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print([DataFolder 'figures/FAU/eleprices'],Format,Resolution);


%% soc

Format='-depsc2';
Resolution=[];

figure('Units','centimeters','Position',[10,7,10,7])
v=linspace(0,24*7,length(S(1).soc));
plot(v,[S(1:3).soc],'LineWidth',1.5)
xlim([0,24*7])
ylim([0.2,1])
xlabel('hour');
ylabel('average state of charge')
legend({'no CS','CS','CS+FCR'},'Orientation','Horizontal','Location','South')
box on
grid on
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print([DataFolder 'figures/FAU/soc'],Format,Resolution);


%% requirements

plotta=1;
testFCRreq(Pmat{3},S(3).sc,S(3).socv,plotta);



return


%% sensitivity over 1 day

ContrVec=[0 1 5 10 20];
mVec=2500:500:5000;

[P,R]=generateplotline3('Munich_clustered',[],'Operations.maxsoc',0.9,'gridday',135,'FCR',FCR,'FCR.contracted',ContrVec,'m',mVec);

FCRfails=zeros(length(ContrVec),length(mVec));

for j=1:size(R,2)
    
    for k=1:size(R,1)
        
        FCRres=testFCR(P{j},R(j));
        FCRfails(k,j)=FCRres.FailMinutes;
        
    end
end


%% plots for 1 day

linevect={'x-','s-','v-','^-','d-','>--','<--','o--','x--'};
% Format='-depsc2';
% Resolution=[];
Format='-dpng';
Resolution='-r300';

figure('Units','centimeters','Position',[10,7,10,7])
hold on
for i=1:size(FCRfails,2)
    plot(FCRfails(:,i),linevect{i})
end
xticks(1:size(FCRfails,1))
xticklabels(num2str(ContrVec'))
xlabel('FCR contracted (MW)')
ylabel('missed 1-minute intervals')
legend({'m=3000','m=4000','m=5000'})
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print([DataFolder 'figures/FAU/FCR1'],Format,Resolution);


%% power plots

plotta(R(1,1),'power','FAU','_m3k_0MW')
plotta(R(5,1),'power','FAU','_m3k_20MW')


%% status

% Res=R(5,1);
Res=Res1;
tsim=Res.Params.tsim;
NumIdle=(sum(Res.Internals.sc(1:tsim,:),2));
x=linspace(0,24,tsim);
xt=0:4:24;

linevect={'x-','s-','v-','^-','d-','>--','<--','o--','x--'};
% Format='-depsc2';
% Resolution=[];
Format='-dpng';
Resolution='-r300';
figure('Units','centimeters','Position',[10,7,10,7])
plot(x,NumIdle)
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('idle vehicles')
legend({'m=3000','m=4000','m=5000'})
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
% print([DataFolder 'figures/FAU/status'],Format,Resolution);



return

%%

figure
plot(sum(Res2.Sim.ef,2))

figure
plot(sum(Res2.Sim.e+Res2.Sim.ef,2))

DataFolder=setDataFolder();


% NOTE: temporary!!!
load([DataFolder 'grid/frequency_30s2015.mat'],'f30s2015');
f=average2(f30s2015,2);
maxf=50.2;
minf=49.8;
af=10000;    % FCR rate per time step (normalized)
FCR=af*min(1,max(-1,(1-(f-minf)/(maxf-minf)*2))); % needed FCR
FCRe=FCR(1:1440);
FCRreal=sum(Res2.Sim.ef,2);

figure
hold on
plot(FCRe)
plot(FCRreal)
plot(FCRreal-FCRe)



% need to check if effective

% plot ideal line
