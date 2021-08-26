%% plot a map of the scenario with charging stations highlighted 

function plotmapcs(scenario)

[~,~,clusters,clusterIDs,chargingStations,C]=getscenario(scenario);

cspos=chargingStations(:,1);

figure
hold on
scatter(C(:,1),C(:,2))
scatter(C(cspos,1),C(cspos,2),'x')
text(C(cspos,1),C(cspos,2),num2str(cspos))

end