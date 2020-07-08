

%% plot of map

DataFolder=setDataFolder();

% load([DataFolder 'scenarios/NewYork/NYC2016/NYGPS200nodes.mat'],'idx','C','SUMS');
load('data/scenarios/NYC2016.mat','C');
load([DataFolder 'scenarios/NewYork/NYC2016/ManhattanMainBoundaries.mat'],'ManhattanMainBoundariesSimplified');
% load([DataFolder 'trips/NYC2016_Jan13-Mar16_10days.mat'])
limiti=[[2;4],[3;6]];
Boundaries=convertToCartesian(ManhattanMainBoundariesSimplified,ManhattanMainBoundariesSimplified,40);
figure
hold on
plot(Boundaries(:,1),Boundaries(:,2))
scatter(C(:,1),C(:,2),'.')
% scatter(A{1}(:,1),A{1}(:,2))
xlabel('km')
ylabel('km')
line([limiti(1,1);limiti(2,1);limiti(2,1);limiti(1,1);limiti(1,1)],[limiti(1,2);limiti(1,2);limiti(2,2);limiti(2,2);limiti(1,2)],'Color','k','LineStyle','--')
axis equal tight
box on
grid on
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
camroll(-90)
print([DataFolder 'figures/Energy/mapNYC200'],'-depsc')