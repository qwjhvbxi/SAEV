

function [dkemd,dkod,dktrip,fk]=generatetripdata3(fo,fd,dk,T,mtsim)

addpath('functions/emd/');
n=size(T,1);

% number of departures and arrival at each station (from original distribution)
fko=zeros(mtsim,n);
fkd=zeros(mtsim,n);
dkod=zeros(mtsim,1);

% length of each interval
tlen=round(1440/mtsim);

for tind=1:mtsim
    
    hind=(tind-1)*tlen+1:tind*tlen;
    fko(tind,:)=sum(fo(hind,:));
    fkd(tind,:)=sum(fd(hind,:));
    dkod(tind)=sum(dk(hind));
    
end

% number of requests per step
fk=sum(fko,2);

% number of vehicle-steps traveled at each macrostep for passengers
e=zeros(mtsim,1);
Flow=cell(mtsim,1);

for tindex=1:mtsim

    % number of vehicle-steps traveled at each macrostep for passengers
    % dkod(tindex)=sum(sum(a(:,:,tindex).*P.t));

    [e(tindex),Flow{tindex}]=emd_mex(fkd(tindex,:),fko(tindex,:),T);

end

dkemd=e.*fk;

% vehicle steps traveled
dktrip=dkod+dkemd;
dktrip(isnan(dktrip))=0;

% dktrip=min(dktrip,P.m*P.macrosteps);

dktrip=repmat(dktrip,2,1);

end


