%% generateClusters(scenario,K[,J,kmed])
% generate a new scenario file with K clusters and J charging stations
% found with k-medoids (kmed=1, default) or k-means (kmed=0). The new
% scenario is named 'scenario_K-J'. If only K is provided, charging
% stations are placed in the cluster centers.

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
    [IDX,~,~,~,nodeIdClusters]=kmedoids(tempC,K);
    
    if K==J
        nodeIdCs=nodeIdClusters;
    else
        % generate charging stations
        [~,~,~,~,nodeIdCs] =kmedoids(tempC,J);
        %     [~,CS2,~,D,nodeIdCs0] =kmedoids(CS,J);
        %     nodeIdCs=nodeIdClusters(nodeIdCs0)    
    end
    
else % k-means
    
    % generate clusters
    [IDX,~,~,D] =kmeans(tempC,K);
    [~,nodeIdClusters]=min(D); % find closest node to station
    
    if K==J
        nodeIdCs=nodeIdClusters;
    else
        % generate charging stations
        [~,~,~,D] =kmeans(tempC,J);
        [~,nodeIdCs]=min(D); % find closest node to station
    end
    
end

% associate to real node number
chargingStations=realT(nodeIdCs);
clusterIDs=realT(nodeIdClusters);
Clusters=zeros(length(T),1);
for i=1:length(IDX)
    Clusters(realT(i))=IDX(i);
end

save([DataFolder 'scenarios/' scenario '_' num2str(K) '-' num2str(J)],'Clusters','clusterIDs','chargingStations')

