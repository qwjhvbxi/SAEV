% P=cpar('NYC2016-20clusters');
addpath functions utilities plots
DataFolder=setDataFolder();
P=cpar('Munich_clustered');
P.Operations.maxsoc=0.9;

FCR.filename='frequency_test';
FCR.limits=[49.8,50.2]; % frequency band
FCR.contracted=5;      % MW
FCR.slowchargeratio=0.8;
FCR.fastchargesoc=0.5;

%% test frequency file

if 0
    testfilename='frequency_test';
    f=repelem(repmat([eye(96)*(FCR.limits(1)-50) eye(96)*(FCR.limits(2)-50)],1,2),15,1)+50;
    save([DataFolder 'grid/' testfilename],'f'); % resolution must be at least 1 minute
end

%% launch test

P.FCR=FCR;
P.m=3000;
P.gridday=17*4+96; % positive from 18:00
Res1=generalC(P,-1,2);

figure
plot(sum(Res1.Sim.ef,2))

[FailMinutes,DeltaPower]=testFCR(P,Res1)

