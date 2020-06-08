%% [A,Atimes,ASortInd,AbuckC, ...
% ODistToNode,ONodeID,DDistToNode,DNodeID,RawDistance]=generateGPStrips(P)
% 
% load trip files and initialize trip variables.
% P contains 'tripfile' and 'scenarioid'. In case of GPS coordinates, P
% contains also the station coordinates 'coords' and (optionally)
% 'stationradius'

function [A,Atimes,ASortInd,AbuckC,ODistToNode,ONodeID,DDistToNode,DNodeID,RawDistance]=generateGPStrips(P)

if ~isfield(P,'tau')
    P.tau=1;
end
if ~isfield(P,'scenarioid') && isfield(P,'tripday') % legacy
    P.scenarioid=P.tripday;
end

DataFolder=setDataFolder();

% retrieve trip matrix
tripFileLocation=[DataFolder 'trips/' P.tripfile '.mat'];
if exist(tripFileLocation,'file') 
    load(tripFileLocation,'A','Atimes');
else
    error('File ''%s'' does not exist in ''%s''.',[P.tripfile '.mat'],[DataFolder 'trips/'])
end
if iscell(A) && P.scenarioid>0 % file with multiple scenarios?
    A1=double(A{P.scenarioid});
    Atimes1=double(Atimes{P.scenarioid});
else
    A1=double(A);
    Atimes1=double(Atimes);
end

% remove 0s
Atimes1(Atimes1(:,1)==0)=1440;

% fast method for trip generation
[Atimes,ASortInd]=sortrows(Atimes1,1); % trips sorted by departure time
A=A1(ASortInd,:); % trips sorted by departure time
Abuck=histc(Atimes(:,1),1:P.tau:1440*P.tau+1); % number of trips in each minute
AbuckC=[0;cumsum(Abuck)]; % total number of trips up to each minute

% GPS coordinates mode
if size(A,2)==4 
    
    % walking speed
    WalkingSpeed=4; % km/h
    
    % check P.coords
    if ~isfield(P,'coords') || length(P.coords)<2
        try 
            load(['data/scenarios/' P.scenario '.mat'],'T','C');
            P.coords=C;
        catch
            error('Must specify P.coords!')
        end
    end

    % raw distance
    RawDistance=sqrt((A(:,1)-A(:,3)).^2+(A(:,2)-A(:,4)).^2);
    
    % distance of each trip from stations
    DistancesToNodesFile=[DataFolder 'temp/' P.tripfile '_' num2str(length(P.coords)) 'DistancesToNodes_Day' num2str(P.scenarioid) '.mat'];
    if exist(DistancesToNodesFile,'file')
        load(DistancesToNodesFile,'ODistToNode','ONodeID','DDistToNode','DNodeID');
    else
        [ODistToNode,ONodeID]=coordsToNodes(A(:,1:2),P.coords,10);
        [DDistToNode,DNodeID]=coordsToNodes(A(:,3:4),P.coords,10);
        save(DistancesToNodesFile,'ODistToNode','ONodeID','DDistToNode','DNodeID');
    end

    % equivalent station matrix (closest station to each request GPS origin)
    A=[ONodeID(:,1),DNodeID(:,1)];

    if isfield(P,'stationradius') 
        
        % only select stations within a radius from each request
        ONodeID(ODistToNode>P.stationradius)=NaN;
        DNodeID(DDistToNode>P.stationradius)=NaN;

        % trim matrices to maximum number of stations within the radius
        Nmax=max(max(sum(~isnan(ONodeID),2)),max(sum(~isnan(DNodeID),2)));
    
    else
        
        % only consider closest station
        Nmax=1;
        
    end
    
    ONodeID(:,Nmax+1:end)=[];
    DNodeID(:,Nmax+1:end)=[];
    
    ONodeID(:,1)=A(:,1);
    DNodeID(:,1)=A(:,2);
    
    % transform distance to time
    ODistToNode=ceil(ODistToNode/WalkingSpeed*60);
    DDistToNode=ceil(DDistToNode/WalkingSpeed*60);
    
else
    ODistToNode=[];
    ONodeID=[];
    DDistToNode=[];
    DNodeID=[];
    RawDistance=[];
end