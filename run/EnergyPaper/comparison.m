%% comparison between optimal and simplified

addpath plots

P1=cpar('NYC2016-small','opti');
Res1=generalC(P1,1,2)

P2=cpar('NYC2016-small','simplified');
P2.Operations.maxwait=Inf;
Res2=generalC(P2,1,2)

[P,R]=generateplotline3('NYC2016',[],'Operations.maxwait',[10 Inf]);




% %% debug simulation
% 
% P1=cpar('NYC2016-small2','simplified');
% P1.Operations.maxwait=Inf;
% Res1=generalC(P1,2,2)
% P=cpar('NYC2016-small2','opti');
% Res2=generalC(P,2,2)
% 
% return