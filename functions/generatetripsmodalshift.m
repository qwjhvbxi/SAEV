%% generate trips depending on modal shift
% With dependence on trip distance
% OriginalModes={'bicycle','car','foot','public_transport'}
% OriginalModes={'bicycle','car','foot','public_transport','car passenger'}
% ModeShift

function generatetripsmodalshift(oldTripName,newTripName,tripDay,ModeShift,D)

addpath functions utilities
DataFolder=getdatafolder();

originalFile=[DataFolder 'trips/' oldTripName '/d' num2str(tripDay)];
load(originalFile,'A','Atimes');

if nargin==5
    res=0.5; % resolution
    tripDistancesKm=D(sub2ind(size(D),A(:,1),A(:,2)))/1000;
    distIndex=min(size(ModeShift,2),max(1,floor(tripDistancesKm/res)));
    IND=sub2ind(size(ModeShift),A(:,3),distIndex);
    
    % select trips depending on mode and distance
    ToSelect=rand(length(IND),1)<ModeShift(IND);
else
    % no distance dependence
    ModeShift=ModeShift(:,1);
    ToSelect=rand(length(A),1)<ModeShift(A(:,3));
end

A=A(ToSelect,:);
Atimes=Atimes(ToSelect,:);

[A,Atimes,~]=cleantripdata(A,Atimes,1);

newFolder=[DataFolder 'trips/' newTripName];
if ~exist(newFolder,'dir')
    mkdir(newFolder)
end

save([newFolder '/d' num2str(tripDay)],'A','Atimes');

end


%     % remove trips that belong to non existent nodes
%     ToRemove=logical(sum(A(:,1)==TzerosID',2)+sum(A(:,2)==TzerosID',2));
%     A(ToRemove,:)=[];
%     Atimes(ToRemove,:)=[];