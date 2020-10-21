% P=cpar('NYC2016-20clusters');
addpath functions utilities plots
DataFolder=setDataFolder();
P=cpar('Munich_clustered_week');

P.Operations.maxsoc=0.8;

FCR.limits=[49.8,50.2]; % frequency band
FCR.contracted=5;      % MW
FCR.slowchargeratio=0.8;
FCR.fastchargesoc=0.4;
FCR.aggregatechargeratio=0.5;
P.FCR=FCR;


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
P.m=5000;
P.gridday=9+96; % positive from 2:00
% P.gridday=17*4+96; % positive from 18:00
% P.gridday=17*4; % negative from 18:00
% P.gridday=16*4+96; % positive from 17:00
Res1=generalC(P,1,2);

figure
plot(sum(Res1.Sim.ef,2))

figure
hold on
plot(sum(Res1.Sim.e+double(Res1.Sim.ef),2),'k-')
plot(sum(Res1.Sim.e,2),'--')
plot(sum(Res1.Sim.ef,2),'r:')


FCRres=testFCR(P,Res1)

