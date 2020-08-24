


function generateClusters(scenario,K,Plots)

addpath functions
DataFolder=setDataFolder();

if nargin==2

    load([DataFolder 'scenarios/' scenario],'T','C')

    % find clusters
    [Clusters,CS]=kmeans(C,K);

    % find closest nodes to cluster centroids
    distances=(C(:,1)-CS(:,1)').^2+(C(:,2)-CS(:,2)').^2;
    [~,chargingStations]=min(distances);
    chargingStations=chargingStations';

    save([DataFolder 'scenarios/' scenario '-' num2str(K) 'clusters'],'T','C','Clusters','chargingStations')

else
    
    if Plots
    
        load([DataFolder 'scenarios/' scenario '-' num2str(K) 'clusters'],'T','C','Clusters','chargingStations')

        c=2
        
        figure
        hold on
        scatter(C(:,1),C(:,2))
%         scatter(CS(:,1),CS(:,2))
        scatter(C(chargingStations,1),C(chargingStations,2),'filled')
        scatter(C(Clusters==c,1),C(Clusters==c,2),'s')
        axis equal tight
    
    end
    
end


