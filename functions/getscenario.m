%% [T,D,clusters,chargingStations]=GETSCENARIO(scenario)
% Retrieve scenario files.
% 
% See also: mainsim

function [T,D,clusters,clusterIDs,chargingStations]=getscenario(scenario)

DataFolder=getdatafolder();

%% load distance matrix

load([DataFolder 'scenarios/' scenario '.mat'],'T','C','D','Clusters','clusterIDs','chargingStations');


%% setup clustering

if exist('Clusters','var')
    clusters=Clusters;
end

if ~exist('T','var')
    k=strfind(scenario,'_');
    scenarioRoot=scenario(1:k(end)-1);
    load([DataFolder 'scenarios/' scenarioRoot '.mat'],'T','C','D');
end

% number of nodes 
if ~isstruct(T)
    n=size(T,1);
else
    n=size(T(1).traveltime,1);
end

if ~exist('Clusters','var')
    chargingStations=(1:n)';
    clusters=(1:n)';
end

if size(chargingStations,2)==1
    chargingStations=[chargingStations , Inf(length(chargingStations),1)];
end

if ~exist('clusterIDs','var')
    clusterIDs=chargingStations(:,1);
end


%% distance

if ~exist('D','var') 
%     if exist('C','var')
%         D=[]; % TODO: calculate approximate distance from coordinates
%     else
    [thisT]=gettraveltimenow(T,0);
    D=30*(thisT/60)*1000; % distance (meters) with average speed of 30 km/h
end





