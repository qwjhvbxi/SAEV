%% generateClusters(scenario,K,Plots)
% generate a new scenario file with K charging stations/clusters found with
% k-means. The new scenario is named 'scenario-Kclusters'
% If Plots==true, it plots the corresponding scenario

function generateClusters(scenario,K,Plots)

addpath functions
DataFolder=getdatafolder();

if nargin<3 || ~Plots

    load([DataFolder 'scenarios/' scenario],'T','C')
    
    T=gettraveltimenow(T,0);
    
    % remove nan values if they exist
    nanT=isnan(T(:,1));
    realT=find(1-nanT);
    tempC=C;
    tempC(nanT,:)=[];

    % generate charging stations
    [IDX,CS,~,D] =kmeans(tempC,K);

    % find closest node to station
    [~,nodeID]=min(D);

    % associate to real node number
    chargingStations=realT(nodeID);
    Clusters=zeros(length(T),1);
    for i=1:length(IDX)
        Clusters(realT(i))=IDX(i);
    end

    T(isnan(T))=0;

    save([DataFolder 'scenarios/' scenario '_' num2str(K)],'Clusters','chargingStations')

else
    
    load([DataFolder 'scenarios/' scenario],'C');
    load([DataFolder 'scenarios/' scenario '_' num2str(K)],'Clusters','chargingStations')

    c=2

    figure
    hold on
    scatter(C(:,1),C(:,2))
    scatter(C(chargingStations,1),C(chargingStations,2),'filled')
    scatter(C(Clusters==c,1),C(Clusters==c,2),'s')
    axis equal tight
    
end


