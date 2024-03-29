%% [T,D,clusters,clusterIDs,chargingStations]=GETSCENARIO(scenario)
% Retrieve scenario files.
% 
% See also: mainsim

function [T,D,clusters,clusterIDs,chargingStations,C]=getscenario(scenario)

DataFolder=getdatafolder();

%% load distance matrix

load([DataFolder 'scenarios/' scenario '.mat'],'T','C','D','Clusters','clusterIDs','chargingStations');

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


%% setup charging stations

if ~exist('chargingStations','var')
    chargingStations=(1:n)';
end

if size(chargingStations,2)==1
    chargingStations=[chargingStations , Inf(length(chargingStations),1)];
end


%% setup clustering

if exist('Clusters','var')
    clusters=Clusters;
    if ~exist('clusterIDs','var') % legacy
        clusterIDs=chargingStations(:,1);
    end
else
    clusters=(1:n)';
    clusterIDs=(1:n)';
end


%% distance

if ~exist('D','var') 
%     if exist('C','var')
%         D=[]; % TODO: calculate approximate distance from coordinates
%     else
    [thisT]=gettraveltimenow(T,0);
    avgspeedkmh=30; % average speed of 30 km/h
    D=avgspeedkmh*(thisT/60); % distance (meters) with average speed 
end





