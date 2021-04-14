%% [A,Atimes,AbuckC,Distances]=gettrips(P)
% P struct contains fields 'tripfolder' and 'tripday' OR 'tripfile'.
% Files should be in folder trips/
%
% See also: main

function [A,Atimes,AbuckC,Distances]=gettrips(P)

%% loading

% set external data folder
DataFolder=getdatafolder();

% determine file name
if isfield(P,'tripfolder') && ~isempty(P.tripfolder)
    if isfield(P,'tripday') && ~isempty(P.tripday)
        
        % tripfolder represents a folder with days inside
        tripFileLocation=[DataFolder 'trips/' P.tripfolder '/d' num2str(P.tripday) '.mat'];
        
        % load files
        [A,Atimes,Distances]=loadFile(tripFileLocation);
        
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
    %tripFileLocation=[DataFolder 'trips/' P.tripfile '.mat'];
    
end


A=[A;A2];
Atimes=[Atimes;(1440+Atimes2)];


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

% load files
if exist(tripFileLocation,'file')
    load(tripFileLocation,'A','Atimes','Distances');
    
    if ~exist('Distances','var')
        Distances=[];
    end
    
else
    %error('File ''%s'' does not exist in ''%s''.',[P.tripfile '.mat'],[DataFolder 'trips/'])
    %     A=NaN;
    %     Atimes=NaN;
    %     Distances=NaN;
    error('File ''%s'' does not exist.',tripFileLocation);
end

end



