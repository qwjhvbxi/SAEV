%% ACartesian=GEO2CARTESIAN(Ageodesic,BoundariesGeodesic,ApproxLatitude)
% Convert geodesic coordinates to approximate cartesian coordinates (km)
% Ageodesic=[Lon,Lat] (Nx2)
% BoundariesGeodesic=[min Lat, min Lon] (1x2)
% ApproxLatitude=Lat (1x1) 

function [ACartesian,BoundariesGeodesic]=geo2cartesian(Ageodesic,BoundariesGeodesic,ApproxLatitude)
    if nargin<3
        ApproxLatitude=mean(Ageodesic(:,2));
    end
    if nargin<2
        BoundariesGeodesic=min(Ageodesic,[],1);
    end
    ACartesian=[(deg2km(Ageodesic(:,1)-min(BoundariesGeodesic(:,1))))*cos(deg2rad(ApproxLatitude)), deg2km(Ageodesic(:,2)-min(BoundariesGeodesic(:,2)))];
end
