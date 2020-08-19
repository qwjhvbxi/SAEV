
if ~exist('P','var')
    
P=cpar('NYC2016');

addpath functions
DataFolder=setDataFolder();
load([DataFolder 'scenarios/' P.scenario],'C')

% find clusters
K=20;
[Clusters,CS]=kmeans(C,K);

% find closest nodes to cluster centroids
distances=(C(:,1)-CS(:,1)').^2+(C(:,2)-CS(:,2)').^2;
[~,chargingStations]=min(distances);

% figure
% hold on
% axis equal
% scatter(C(:,1),C(:,2))
% scatter(CS(:,1),CS(:,2))
% scatter(C(chargingStations,1),C(chargingStations,2),'x')
% scatter(C(Clusters==1,1),C(Clusters==1,2),'s')

P.chargingStations=chargingStations';
P.clusters=Clusters;
P.Operations.maxidle=5; % minutes

end

%% 

P.e=1;

Res=generalC(P,-1,2)


%% check for contraint violations

if 0

    d=Res.Internals.d;
    
    % charging when not in charging stations and with delay=0
    v1=(d>0).*Res.Internals.s2;
    sum(v1(:))

    % charging more than possible

% charging with wrong status

% moving while charging

% delay
% should be not more than 1
dchanges=max(Res.Internals.d(1:end-1,:)-Res.Internals.d(2:end,:));
find(dchanges>1)

v=9972;
z=1:100;
% d s1 s2 s3
full([Res.Internals.d(z,v) Res.Internals.s1(z,v) Res.Internals.s2(z,v) Res.Internals.s3(z,v)])
plot(Res.Internals.d(z,v))

end



