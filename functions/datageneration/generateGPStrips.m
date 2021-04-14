%% generateGPStrips(FileIn,FileOut,C,stationradius,WalkingSpeed)
% create trip files from coordinates
%
% See also generalC

function generateGPStrips(FileIn,FileOut,C,StationRadius,WalkingSpeed)

if nargin<5    
    % walking speed
    WalkingSpeed=4; % km/h
end

if nargin<4
    StationRadius=[];
end

% load coordinate file
load(FileIn,'A','Atimes');

% clean up data
[A,Atimes,~]=cleantripdata(A,Atimes);
A=double(A);
Atimes=double(Atimes);

% GPS coordinates mode
if size(A,2)==4 

    % raw distance
    Distances.RawDistance=sqrt((A(:,1)-A(:,3)).^2+(A(:,2)-A(:,4)).^2);
    
    [ODistToNode,ONodeID]=computedistancetonodes(A(:,1:2),C,10);
    [DDistToNode,DNodeID]=computedistancetonodes(A(:,3:4),C,10);

    % equivalent station matrix (closest station to each request GPS origin)
    A=[ONodeID(:,1),DNodeID(:,1)];

    if ~isempty(StationRadius) 
        
        % only select stations within a radius from each request
        ONodeID(ODistToNode>StationRadius)=NaN;
        DNodeID(DDistToNode>StationRadius)=NaN;

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
    
    Distances.ONodeID=ONodeID;
    Distances.DNodeID=DNodeID;
    
    
    % transform distance to time
    Distances.ODistToNode=ceil(ODistToNode/WalkingSpeed*60);
    Distances.DDistToNode=ceil(DDistToNode/WalkingSpeed*60);
    
else
    Distances=[];
end

% remove trips inside same node
SameNode=logical(A(:,1)==A(:,2));
A(SameNode,:)=[];
Atimes(SameNode,:)=[];
if ~isempty(Distances) && sum(SameNode)>0
    Distances.ODistToNode(SameNode,:)=[];
    Distances.DDistToNode(SameNode,:)=[];
    Distances.ONodeID(SameNode,:)=[];
    Distances.DNodeID(SameNode,:)=[];
    Distances.RawDistance(SameNode,:)=[];
end

% save trip file
save(FileOut,'A','Atimes','Distances');