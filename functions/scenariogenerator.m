
%% scenario characteristics
scenarioname='Tokyo5x5'; % small map (5x5 km)
tripsperday=200;
n=10;
limiti=[20 25 22 27];
speedkmh=30;

%% generate scenario data
[OD,TPH,T]=generateTokyoScenario(scenarioname,limiti,n,speedkmh,0);

%% generate passenger requests
% retrieve distribution from seed or create new if non existant
% each seed is associated with two distributions. The first one is used by
% the optimization. In case of perfect knowledge it is also used in the
% simulation, otherwise the second is used in the simulation.
namearrivals=['data/scenarios/trips/' scenarioname '-' num2str(tripsperday) '-' num2str(1) '.mat'];
if exist(namearrivals,'file')
    load(namearrivals,'c1','c2')
else
    c1=generatearrivals2(OD,TPH,n,tripsperday);
    c2=generatearrivals2(OD,TPH,n,tripsperday);
    save(namearrivals,'c1','c2');
end