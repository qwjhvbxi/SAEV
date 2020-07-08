%% convert geodesic coordinates to approximate cartesian coordinates (km)
% Ageodesic=[Lon,Lat] (Nx2)
% BoundariesGeodesic=[min Lat, min Lon] (1x2)
% ApproxLatitude=Lat (1x1) 

function ACartesian=convertToCartesian(Ageodesic,BoundariesGeodesic,ApproxLatitude)
    ACartesian=[(deg2km(Ageodesic(:,1)-min(BoundariesGeodesic(:,1))))*cos(deg2rad(ApproxLatitude)), deg2km(Ageodesic(:,2)-min(BoundariesGeodesic(:,2)))];
end
