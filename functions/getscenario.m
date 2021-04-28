function [T,clusters,chargingStations]=getscenario(P)

DataFolder=getdatafolder();

%% load distance matrix

load([DataFolder 'scenarios/' P.scenario '.mat'],'T','Clusters','chargingStations');


%% setup clustering

if exist('Clusters','var')
    clusters=Clusters;
end

if ~exist('T','var')
    k=strfind(P.scenario,'_');
    scenarioRoot=P.scenario(1:k(end)-1);
    load([DataFolder 'scenarios/' scenarioRoot '.mat'],'T');
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

