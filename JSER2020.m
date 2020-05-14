%% simulations for JSER paper



%% generate Tokyo scenario 

DataFolder=setDataFolder();

load([DataFolder 'input_files/Tokyo/sources/tripsCell.mat'],'Trips');
load([DataFolder 'input_files/Tokyo/sources/aev18data.mat'],'areas','chargingstnodes','gridcoord','tph');

% each trip in Trips is generated according to a gaussian centered in the
% hour. Weights are between 16 and 210. I divide by 15 and round, so I get
% about 188,000 trips

%% create A, Atimes

ScaleFactor=15;
A=[];
Atimes=[];
for h=1:24
    
    TripWeights=max(1,round(Trips{h}(:,1)./ScaleFactor));
    HourlyTrips=repelem(Trips{h}(:,2:3),TripWeights,1);
    RandMinutes=round((h-1)*60+30+randn(length(HourlyTrips),1)*30);
    A=[A;HourlyTrips];
    Atimes=[Atimes;[RandMinutes RandMinutes]];
    
end

% bring trips out of either side of the day on the other side
Atimes(Atimes>1440)=Atimes(Atimes>1440)-1440;
Atimes(Atimes<1)=Atimes(Atimes<1)+1440;

% remove trips starting and ending in same node
SameNode=A(:,1)==A(:,2);
A(SameNode,:)=[];
Atimes(SameNode,:)=[];

% save([DataFolder 'trips/Tokyo2008_1day_189k.mat'],'A','Atimes'); % with same node trips
save([DataFolder 'trips/Tokyo2008_1day_162k.mat'],'A','Atimes');


%% create scenario with less stations

limits=[10 30 15 35];
selection=logical((gridcoord(:,1)>limits(1)).*(gridcoord(:,1)<limits(2)).*(gridcoord(:,2)>limits(3)).*(gridcoord(:,2)<limits(4)));
AllowedStations=find(selection);

% create vector to assign new relative indeces to nodes
NewIndeces=zeros(514,1);
NewIndeces(AllowedStations)=1:length(AllowedStations);

% figure
% scatter(gridcoord(:,1),gridcoord(:,2))
% figure
% scatter(gridcoord(selection,1),gridcoord(selection,2))

TripSelection=logical(ismember(A(:,1),AllowedStations).*ismember(A(:,2),AllowedStations));
A=NewIndeces(A(TripSelection,:));
Atimes=Atimes(TripSelection,:);
save([DataFolder 'trips/Tokyo2008_1day_48k.mat'],'A','Atimes');


%% create C, T (from Euclidean distance)

Eta=1.5;
T=max(1,round(Eta*sqrt((gridcoord(:,1)-gridcoord(:,1)').^2+(gridcoord(:,2)-gridcoord(:,2)').^2)));

C=gridcoord;
save('data/scenarios/Tokyo514.mat','C','T');

C=gridcoord(selection,:);
T=T(selection,selection);
save('data/scenarios/Tokyo189.mat','C','T')








%% simulations