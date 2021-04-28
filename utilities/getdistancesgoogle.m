%% [d,resvec]=GETDISTANCESGOOGLE(apikey,lat1,lon1,lat2,lon2,traveldate)
% Querying Google Maps Distance Matrix API
% traveldate: dd/MM/yyyy HH:mm

function [d,resvec]=getdistancesgoogle(apikey,lat1,lon1,lat2,lon2,traveldate)

addpath('utilities/jsonlab');

% input formatting
originvec=reshape([lat1,lon1]',[1,length(lat1)*2])';
destinvec=reshape([lat2,lon2]',[1,length(lat2)*2])';
originstr=sprintf('%3.4f,%3.4f|',originvec);
destinstr=sprintf('%3.4f,%3.4f|',destinvec);
origintime=num2str(posixtime(datetime(traveldate,'inputFormat','dd/MM/yyyy HH:mm')));

% url1 = ['https://maps.google.co.jp/maps/api/geocode/json?address=' , urlencode([Codesraw.prefectures{i}{1} , Codesraw.wards{i}{1} , addr1{k}]) , '&sensor=false&language=ja&region=JP&key=' apikey];
url1 = ['https://maps.google.co.jp/maps/api/distancematrix/json?origins=',originstr,'&destinations=',destinstr,'&departure_time=',origintime,'&mode=driving&sensor=false&language=ja&region=JP&key=' apikey];
str1 = urlread(url1,'Charset','UTF-8');
resvec = loadjson(str1);

% if result found
if strcmp(resvec.status,'OK')
    d=resvec.rows{1}.elements{1}.duration.value;
else
    d=NaN;
end

end







