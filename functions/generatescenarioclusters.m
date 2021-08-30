%% generateClusters(scenario,K,J[,kmed])
% generate a new scenario file with J charging stations and K clusters found with
% k-medoids (kmed=1) or k-means. The new scenario is named 'scenario_J-K'

function generatescenarioclusters(scenario,K,J,kmed)

if nargin<3
    J=K;
end
if nargin<4
    kmed=true;
end

addpath functions
DataFolder=getdatafolder();

load([DataFolder 'scenarios/' scenario],'T','C')

T=gettraveltimenow(T,0);

% remove nan values if they exist
nanT=isnan(T(:,1));
realT=find(1-nanT);
tempC=C;
tempC(nanT,:)=[];

if kmed

    % generate clusters
    [IDX,~,~,~,nodeIdClusters] =kmedoids(tempC,K);
    
    % generate charging stations
    [~,~,~,~,nodeIdCs] =kmedoids(tempC,J);
    
%     [~,CS2,~,D,nodeIdCs0] =kmedoids(CS,J);
%     nodeIdCs=nodeIdClusters(nodeIdCs0)
    
else
    
    % generate clusters
    [IDX,~,~,D] =kmeans(tempC,K);
    [~,nodeIdClusters]=min(D); % find closest node to station
    
    % generate charging stations
    [~,~,~,D] =kmeans(tempC,J);
    [~,nodeIdCs]=min(D); % find closest node to station
    
end

% associate to real node number
chargingStations=realT(nodeIdCs);
clusterIDs=realT(nodeIdClusters);
Clusters=zeros(length(T),1);
for i=1:length(IDX)
    Clusters(realT(i))=IDX(i);
end

save([DataFolder 'scenarios/' scenario '_' num2str(K)],'Clusters','clusterIDs','chargingStations')

