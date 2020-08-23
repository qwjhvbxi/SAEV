


function generateClusters(scenario,K,Plots)

addpath functions
DataFolder=setDataFolder();
load([DataFolder 'scenarios/' scenario],'T','C')

% find clusters
[Clusters,CS]=kmeans(C,K);

% find closest nodes to cluster centroids
distances=(C(:,1)-CS(:,1)').^2+(C(:,2)-CS(:,2)').^2;
[~,chargingStations]=min(distances);
chargingStations=chargingStations';

save([DataFolder 'scenarios/' scenario '-' num2str(K) 'clusters'],'T','C','Clusters','chargingStations')

if nargin>2 && Plots
    figure
    hold on
    axis equal
    scatter(C(:,1),C(:,2))
    scatter(CS(:,1),CS(:,2))
    scatter(C(chargingStations,1),C(chargingStations,2),'x')
    scatter(C(Clusters==1,1),C(Clusters==1,2),'s')
end


