%% [A,Atimes,AbuckC,Distances]=GETTRIPS(P)
% Get trip files from information in parameter struct. 
% P struct contains fields 'tripfolder' and optionally 'tripday'.
% If 'tripday' is not set, then 'tripfolder' is treated as a single file
% instead of a folder.
% All trip files/folders should be in folder 'data/trips/'.
%
% See also: main

function [A,Atimes,AbuckC,Distances]=gettrips(P)

% set external data folder
DataFolder=getdatafolder();

% determine file name
if isfield(P,'tripfolder') && ~isempty(P.tripfolder)
    if isfield(P,'tripday') && ~isempty(P.tripday)
        
        % tripfolder represents a folder with days inside
        tripFileLocation=[DataFolder 'trips/' P.tripfolder '/d' num2str(P.tripday) '.mat'];
        
        % load files
        [A,Atimes,Distances]=loadFile(tripFileLocation);
        
        % load next day if it exists
        tripFileLocation2=[DataFolder 'trips/' P.tripfolder '/d' num2str(P.tripday+1) '.mat'];
        if exist(tripFileLocation2,'file')
            [A2,Atimes2,~]=loadFile(tripFileLocation2);
        else
            A2=A;
            Atimes2=Atimes;
        end
        
    else
        
        % tripfolder is actually a file with a single day
        tripFileLocation=[DataFolder 'trips/' P.tripfolder '.mat'];
        
        % load files
        [A,Atimes,Distances]=loadFile(tripFileLocation);
        
        A2=A;
        Atimes2=Atimes;
        
    end
else
    
    error('need to specify trip folder!')
    
end

A=double([A;A2]);
Atimes=double([Atimes;(1440+Atimes2)]);


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

end



function [A,Atimes,Distances]=loadFile(tripFileLocation)

if exist(tripFileLocation,'file')
    
    load(tripFileLocation,'A','Atimes','Distances');
    
    if ~exist('Distances','var')
        Distances=[];
    end
    
else
    
    error('File ''%s'' does not exist.',tripFileLocation);
    
end

end



