
function [dkemd,dkod,dktrip,fk]=generatetripdataAlt(fo,fd,dk,T,mtsim)

n=size(T,1);

% length of each interval
tlen=round(1440/mtsim);

dkod=sum(reshape(dk(1:1440),tlen,mtsim))';
dkemd=zeros(mtsim,1);

b=zeros(mtsim,n);

for kt=1:mtsim
    
    ThisInterval=(kt-1)*tlen+1:kt*tlen;
 
%     bmin=0;
%     uv=0;
%     dw=0;
%
%     % calculate imbalance
%     b(kt,:)=uv-dw ...  number of vehicles and passengers waiting at each station
%             +sum(fd(ThisInterval,:)) ...  expected arrivals between now and now+P.ts
%             -sum(fo(ThisInterval,:));   % expected requests between now and now+P.ts+t
%         
%     % identify feeder and receiver stations
%     F=(b(kt,:)-bmin).*(b(kt,:)>=bmin); % feeders
%     R=(-b(kt,:)+bmin).*(b(kt,:)<bmin); % receivers
    
    F=sum(fd(ThisInterval,:));
    R=sum(fo(ThisInterval,:));

    % identify optimal relocation flux
    x=optimalrelocationfluxes(F,R,T);

    % read results
    [Fs,Rs,Vr]=find(x);

    % distance of relocation
    dkemd(kt)=sum(T(sub2ind(size(T),Fs,Rs)).*Vr);


end

dktrip=dkemd+dkod;
dktrip=repmat(dktrip,2,1);

fk=0;



