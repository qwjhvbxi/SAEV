%% plot a map of the scenario with charging stations highlighted 

function plotmapcs(scenario)

[~,~,clusters,clusterIDs,chargingStations,C]=getscenario(scenario);

cspos=chargingStations(:,1);

figure
hold on
if length(clusterIDs)<length(clusters)
    for i=1:length(clusterIDs)
        thispos=(clusters==i);
        scatter(C(thispos,1),C(thispos,2))
    end
else
    scatter(C(:,1),C(:,2))
end
scatter(C(cspos,1),C(cspos,2),'x')
text(C(cspos,1),C(cspos,2),num2str(cspos))
axis equal tight
end