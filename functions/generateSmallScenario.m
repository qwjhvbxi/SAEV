

% generate C coordinates: one day is ok, with optional limits (this file)
% generate T times: need as much data as possible
% generate tripfile: need C and T, optional limits

addpath functions

DataFolder=setDataFolder();

P=cpar('NYC2016');

tripFileLocation=[DataFolder 'trips/' P.tripfile '.mat'];
load(tripFileLocation,'A','Atimes');
Ac=A;
Actimes=Atimes;

setlimits=[0,5;0,5];

[T,C,A2,A2times,tripsubset]=generateScenario(Ac{P.scenarioid},Actimes{P.scenarioid},10,setlimits);


Ac10=Ac{P.scenarioid}(tripsubset,:);
colori=lines(10);
figure
hold on
for i=1:10
    scatter(Ac10(A2(:,1)==i,1),Ac10(A2(:,1)==i,2),1,colori(i,:))
end

scatter(C(:,1),C(:,2),5)

axis equal

save('data/scenarios/NYC2016-small.mat','C','T');


%% select only a part of the trips

% remove trips with same origin/destination
DifferentOD=(A2(:,1)~=A2(:,2));
A3=A2(DifferentOD,:);
A3times=A2times(DifferentOD,:);
Ac10_2=Ac10(DifferentOD,:);


k=300;
M=randperm(length(A3),k);
A=A3(M,:);
Atimes=A3times(M,:);


Ac10_3=Ac10_2(M,:);
figure
hold on
for i=1:10
    scatter(Ac10_3(A(:,1)==i,1),Ac10_3(A(:,1)==i,2),5,colori(i,:))
end
scatter(C(:,1),C(:,2),50)
axis equal

figure
plot(histc(Atimes(:,1),0:10:1440))

save([DataFolder 'trips/NYC2016-small_13Jan.mat'],'A','Atimes');


%% even smaller scenario

load([DataFolder 'trips/NYC2016-small_13Jan.mat'],'A','Atimes');

k=100;
M=randperm(length(A),k);
A=A(M,:);
Atimes=Atimes(M,:);
save([DataFolder 'trips/NYC2016-small2_13Jan.mat'],'A','Atimes');
