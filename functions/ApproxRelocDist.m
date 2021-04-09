%% [dkemd,dkod,dktrip]=ApproxRelocDist(A,Atimes,T,Beta)
% calculate approximate relocation distance needed 

function dkemd=ApproxRelocDist(A,Atimes,T,Beta)

n=size(T,1);
t=double(max(Atimes(:)));
mtsim=round(t/Beta);
t=Beta*mtsim;

fo=zeros(t,n);
fd=zeros(t,n);

for i=1:t
    ThisMinute=logical(Atimes(:,1)==i);
    fo(i,:)=accumarray([A(ThisMinute,1);n],[ones(sum(ThisMinute),1);0]);
    fd(i,:)=accumarray([A(ThisMinute,2);n],[ones(sum(ThisMinute),1);0]);
    
end

dkemd=zeros(mtsim,1);

approx=(n>200);

fprintf('\n Calculating approximate relocation distance\n\n')

for kt=1:mtsim
    
    % progress report
    clc
    fprintf('\n %d/%d\n\n',kt,mtsim)
    
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




