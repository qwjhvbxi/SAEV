
addpath functions

% P=cpar('Munich');
% 
% DataFolder=setDataFolder();
% load([DataFolder 'scenarios/' P.scenario],'C')
% 
% % find clusters
% K=20;
% [Clusters,CS]=kmeans(C,K);
% 
% % find closest nodes to cluster centroids
% distances=(C(:,1)-CS(:,1)').^2+(C(:,2)-CS(:,2)').^2;
% [~,chargingStations]=min(distances);
% 
% % figure
% % hold on
% % axis equal
% % scatter(C(:,1),C(:,2))
% % scatter(CS(:,1),CS(:,2))
% % scatter(C(chargingStations,1),C(chargingStations,2),'x')
% % scatter(C(Clusters==1,1),C(Clusters==1,2),'s')
% 
% P.chargingStations=chargingStations';
% P.clusters=Clusters;
% P.Operations.maxidle=5; % minutes

P=cpar('Munich_clustered');
P.m=3000;
P.e=1;
P.Operations.maxwait=Inf;
P.Tech.chargekw=22;
P.gridday=135;
P.gridfile='Germany_DA_2019';
Res=generalC(P,1,2)

plotta(Res,'power','FAU')