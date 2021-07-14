%% [d,resvec]=GETDISTANCESGOOGLE(apikey,lat1,lon1,lat2,lon2,traveldate)
% Querying Google Maps Distance Matrix API
% traveldate: dd/MM/yyyy HH:mm

function [resvec]=getdistancesgoogle(apikey,lat1,lon1,lat2,lon2,traveldate,mode)

if nargin<7
    mode='driving'; % or transit
end

addpath('utilities/jsonlab');

if exist('loadjson','file')~=2
    warning('Need jsonlab (https://uk.mathworks.com/matlabcentral/fileexchange/33381-jsonlab-a-toolbox-to-encode-decode-json-files).')
    resvec=[];
    return
end

% input formatting
originvec=reshape([lat1,lon1]',[1,length(lat1)*2])';
destinvec=reshape([lat2,lon2]',[1,length(lat2)*2])';
originstr=sprintf('%3.4f,%3.4f|',originvec);
originstr=originstr(1:end-1);
destinstr=sprintf('%3.4f,%3.4f|',destinvec);
destinstr=destinstr(1:end-1);
origintime=num2str(posixtime(datetime(traveldate,'inputFormat','dd/MM/yyyy HH:mm')));

% url1 = ['https://maps.google.co.jp/maps/api/geocode/json?address=' , urlencode([Codesraw.prefectures{i}{1} , Codesraw.wards{i}{1} , addr1{k}]) , '&sensor=false&language=ja&region=JP&key=' apikey];
% url1 = ['https://maps.google.co.jp/maps/api/distancematrix/json?origins=',originstr,'&destinations=',destinstr,'&departure_time=',origintime,'&mode=',mode,'&sensor=false&language=ja&region=JP&key=' apikey];
url1 = ['https://maps.google.co.jp/maps/api/distancematrix/json?origins=',originstr,'&destinations=',destinstr,'&departure_time=',origintime,'&mode=',mode,'&key=' apikey];
str1 = urlread(url1,'Charset','UTF-8');
resvec = loadjson(str1);

end







