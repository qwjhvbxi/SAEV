
if ~exist('P','var')
    
P=cpar('NYC2016');

addpath functions
DataFolder=setDataFolder();
load([DataFolder 'scenarios/' P.scenario],'C')

% find clusters
K=20;
[Clusters,CS]=kmeans(C,K);

% find closest nodes to cluster centroids
distances=(C(:,1)-CS(:,1)').^2+(C(:,2)-CS(:,2)').^2;
[~,chargingStations]=min(distances);

% figure
% hold on
% axis equal
% scatter(C(:,1),C(:,2))
% scatter(CS(:,1),CS(:,2))
% scatter(C(chargingStations,1),C(chargingStations,2),'x')
% scatter(C(Clusters==1,1),C(Clusters==1,2),'s')

P.chargingStations=chargingStations';
P.clusters=Clusters;
P.Operations.maxidle=5; % minutes

end

%% 

P.e=1;
P.tripday=13;

[P1,R1]=generateplotline3('NYC2016',[],'Operations.maxwait',Inf,'tripfolder',["NYC2016-small_150","NYC2016-small_100"],'tripday,gridday',Days);
Res=generalC(P,-1,2)
