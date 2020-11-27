
function [dkemd,dkod,dktrip]=generatetripdataAlt(A,Atimes,T,Beta)

n=size(T,1);
t=double(max(Atimes(:,2)));

% cleanup A
Atimes(Atimes(:,1)==0,1)=1;
Atimes(Atimes(:,2)==0,2)=1;

fo=zeros(t,n);
fd=zeros(t,n);
dk=zeros(t,1);

for i=1:t
    ThisMinute=logical(Atimes(:,1)==i);
    fo(i,:)=accumarray([A(ThisMinute,1);n],[ones(sum(ThisMinute),1);0]);
    fd(i,:)=accumarray([A(ThisMinute,2);n],[ones(sum(ThisMinute),1);0]);
    dk(i)=sum(T(sub2ind(size(T),A(ThisMinute,1),A(ThisMinute,2))));
end

% fo=sparse(fo);
% fd=sparse(fd);

mtsim=round(t/Beta);

dkod=sum(reshape(dk(1:t),Beta,mtsim))';
dkemd=zeros(mtsim,1);

approx=(n>200);

for kt=1:mtsim
    
    ThisInterval=(kt-1)*Beta+1:kt*Beta;
     
    F=sum(fd(ThisInterval,:));
    R=sum(fo(ThisInterval,:));

    % identify optimal relocation flux
    x=optimalrelocationfluxes(F,R,T,60,approx);

    % read results
    [Fs,Rs,Vr]=find(x);

    % distance of relocation
    dkemd(kt)=sum(T(sub2ind(size(T),Fs,Rs)).*Vr);

end

dktrip=dkemd+dkod;



