function [dkemd,dkod,dktrip,fk]=generatetripdata2(OD,T,mtsim)

addpath('emd/');
n=size(T,1);

% % number of departures and arrival at each station (from Poisson)
% fkd=squeeze(sum(a(:,:,:),2))';
% fko=squeeze(sum(a(:,:,:),1))';

% number of departures and arrival at each station (from original distribution)
fko=zeros(mtsim,n);
fkd=zeros(mtsim,n);
dkod=zeros(mtsim,1);


for tind=1:mtsim
    
    hind=ceil(tind/(mtsim/length(OD))); % hour
    fko(tind,:)=accumarray([OD{hind}(:,2);n],[OD{hind}(:,1);0])/(mtsim/length(OD));
    fkd(tind,:)=accumarray([OD{hind}(:,3);n],[OD{hind}(:,1);0])/(mtsim/length(OD));
    dkod(tind)=sum(T(sub2ind([n,n],OD{hind}(:,2),OD{hind}(:,3))).*OD{hind}(:,1)/(mtsim/length(OD)));
    
end

% number of requests per step
fk=sum(fkd,2);

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


