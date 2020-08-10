%% create scenario file

DataFolder=setDataFolder();

% load([DataFolder 'trips/NY_scenario.mat'],'T');
load('../model_aev5/data/scenarios/NYC2018/NY_scenario.mat','T');

namefile=['data/scenarios/NY.mat'];
T=floor(T/60);
save(namefile,'T')


%% create arrivals files

load([DataFolder 'trips/NY_trips_10wed_0103-0307.mat'],'A','Atimes','fd','fo');

k=1;

[c1]=convertAtimes(double(A{k}),ceil(double(Atimes{k})/60),double(max(max(A{k}))),1440);
c1=cat(3,c1,c1(:,:,1:180));
c2=c1;
namearrivals=[DataFolder 'trips_saev/NY-' num2str(k) '.mat'];
save(namearrivals,'c1','c2')