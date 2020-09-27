% P=cpar('NYC2016-20clusters');
addpath functions utilities plots
DataFolder=setDataFolder();
P=cpar('Munich_clustered');
P.Operations.maxsoc=0.9;
P.gridday=135;

FCR.filename='frequency_30s2015';
FCR.limits=[49.8,50.2]; % frequency band
FCR.contracted=10;      % MW
FCR.slowchargeratio=0.8;
FCR.fastchargesoc=0.5;

%% test frequency file

testfilename='frequency_test';
f=repelem(repmat([eye(96)*(FCR.limits(1)-50) eye(96)*(FCR.limits(2)-50)],1,2),15,1)+50;
save([DataFolder 'grid/' testfilename],'f'); % resolution must be at least 1 minute

%% launch test

FCRtest=FCR;
FCRtest.filename='frequency_test';
FCRtest.contracted=10;
Ptest=P;
Ptest.FCR=FCRtest;
Ptest.m=3000;
Ptest.gridday=17*4+96; % positive from 18:00
Res1=generalC(Ptest,-1,2)
figure
plot(sum(Res1.Sim.ef,2))

%% expected FCR

load([DataFolder 'grid/' FCR.filename],'f');
ReshapeFactor=size(f,1)/1440*P.e;
f=average2(f(:,P.gridday),ReshapeFactor);
FCRe=min(1,max(-1,(1-(f-FCR.limits(1))/(FCR.limits(2)-FCR.limits(1))*2))); % needed FCR
FCRe=FCRe(1:1440);

%% launching sims

P1=P;
P2=P;

P1.FCR=FCR;
P2.FCR=[];

Res1=generalC(P1,-1,2)
Res2=generalC(P2,-1,2)

addpath plots
plotta(Res1,'power')
plotta(Res2,'power')

%%

% ContrVec=[0 1 5 10 20];
% mVec=[3000,4000,5000];
% FreqExpect=FCRe*ContrVec*1000;

ContrVec=[0];
mVec=[2500:100:3000];

[P,R]=generateplotline3('Munich_clustered',[],'Operations.maxsoc',0.9,'gridday',135,'FCR',FCR,'FCR.contracted',ContrVec,'m',mVec);

FCRfails=zeros(length(ContrVec),length(mVec));

for j=1:size(R,2)

    S=[R(:,j).Sim];
    Freq=reshape(sum(vertcat(S.ef),2),1440,size(R,1));

    for k=1:size(R,1)
        GridIn=R(k,j).Sim.ef-R(k,j).Sim.e;
        Violations(k,j)=sum(sum(abs(GridIn)>Par.Tech.chargekw));
    end
    
    FCRfails(:,j)=(sum(abs(Freq-FreqExpect)>10^-7));
    
end

%% plots

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
