%% [A,Atimes,ASortInd,AbuckC,Distances]=generateGPStrips(P)
% 
% load trip files and initialize trip variables.
% P contains 'tripfile' and 'scenarioid'. In case of GPS coordinates, P
% contains also the station coordinates 'coords' and (optionally)
% 'stationradius'
%
% See also generalC

function [A,Atimes,AbuckC,Distances]=generateGPStrips(P)

% set external data folder
DataFolder=setDataFolder();

% determine file name
if isfield(P,'tripfolder')
    tripFileLocation=[DataFolder 'trips/' P.tripfolder '/d' num2str(P.tripday) '.mat'];
    TripName=P.tripfolder;
else
    tripFileLocation=[DataFolder 'trips/' P.tripfile '.mat'];
    TripName=P.tripfile;
end

% load files
if exist(tripFileLocation,'file')
    load(tripFileLocation,'A','Atimes','Cleaned');
else
    %error('File ''%s'' does not exist in ''%s''.',[P.tripfile '.mat'],[DataFolder 'trips/'])
    error('File ''%s'' does not exist.',tripFileLocation);
end

if exist('Cleaned','var')
    
    % clean up data
    [A,Atimes,~]=cleanData(A,Atimes);
    Cleaned=true;
    save(tripFileLocation,'A','Atimes','Cleaned');
    
end 

A=double(A);
Atimes=double(Atimes);


% GPS coordinates mode
if size(A,2)==4 
    
    % walking speed
    WalkingSpeed=4; % km/h
    
    % check P.coords
    if ~isfield(P,'coords') || length(P.coords)<2
        try
            load(['data/scenarios/' P.scenario '.mat'],'C');
            P.coords=C;
        catch
            error('Must specify P.coords!')
        end
    end

    % raw distance
    RawDistance=sqrt((A(:,1)-A(:,3)).^2+(A(:,2)-A(:,4)).^2);
    
    % distance of each trip from stations
    DistancesToNodesFile=[DataFolder 'temp/' TripName '_' num2str(length(P.coords)) 'DistancesToNodes_Day' num2str(P.tripday) '.mat'];
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

% remove trips inside same node
SameNode=logical(A(:,1)==A(:,2));
A(SameNode,:)=[];
Atimes(SameNode,:)=[];

% fast method for trip generation
Abuck=histc(Atimes(:,1),1:1441); % number of trips in each minute
AbuckC=[0;cumsum(Abuck)]; % total number of trips up to each minute

% NOTE: need to fix this problem with trips in same node
Distances.ODistToNode=ODistToNode;
Distances.ONodeID=ONodeID;
Distances.DDistToNode=DDistToNode;
Distances.DNodeID=DNodeID;
Distances.RawDistance=RawDistance;