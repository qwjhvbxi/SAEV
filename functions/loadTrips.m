%% [A,Atimes,AbuckC,Distances]=loadTrips(P)
% P struct contains fields 'tripfolder' and 'tripday' OR 'tripfile'.
% Files should be in folder trips/
%
% See also: generalC

function [A,Atimes,AbuckC,Distances]=loadTrips(P)

%% loading

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
    load(tripFileLocation,'A','Atimes','Distances');
else
    %error('File ''%s'' does not exist in ''%s''.',[P.tripfile '.mat'],[DataFolder 'trips/'])
    error('File ''%s'' does not exist.',tripFileLocation);
end


%% checks

% checks: nodes
if sum(floor(A(:))~=A(:))>0
    error('nodes should be integers')
end
    
% checks: ordered
if ~issorted(Atimes(:,1)) 
    error('file not sorted')
end

% checks: same node
if sum(A(:,1)==A(:,2))>0
    error('trips to same node')
end

% checks: impossible trips
if sum(Atimes(:)==0)>0
    error('trips should have non-zero starting time')
end


%% utilities

% fast method for trip generation
Abuck=histc(Atimes(:,1),1:1441); % number of trips in each minute
AbuckC=[0;cumsum(Abuck)]; % total number of trips up to each minute

if ~exist('Distances','var')
    Distances=[];
end

end