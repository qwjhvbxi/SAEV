%% [T,clusters,chargingStations]=GETSCENARIO(scenario)
% Retrieve scenario files.
% 
% See also: mainsim

function [T,clusters,chargingStations,D]=getscenario(scenario)

DataFolder=getdatafolder();

%% load distance matrix

load([DataFolder 'scenarios/' scenario '.mat'],'T','C','D','Clusters','chargingStations');


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


%% distance

if ~exist('D','var') 
    if exist('C','var')
        D=[]; % TODO: calculate approximate distance from coordinates
    else
        D=[];
    end
end




