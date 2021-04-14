%% dkemd=COMPUTEEMD(A,Atimes,T,Beta)
% Compute travel time for relocation during time horizon
% 
% See also: main

function dkemd=computeemd(A,Atimes,T,Beta)

fprintf('\n Calculating approximate relocation distance\n\n')

n=size(T,1);
t=double(max(Atimes(:)));
mtsim=round(t/Beta);
approx=(n>100);

fo=zeros(mtsim,n);
fd=zeros(mtsim,n);

for i=1:mtsim
    thisStep=logical((Atimes(:,1)>=(i-1)*Beta+1).*(Atimes(:,1)<=i*Beta));
    fo(i,:)=accumarray([A(thisStep,1);n],[ones(sum(thisStep),1);0]);
    fd(i,:)=accumarray([A(thisStep,2);n],[ones(sum(thisStep),1);0]);
end

dkemd=zeros(mtsim,1);

for kt=1:mtsim
    
    % progress report
    clc
    fprintf('\n %d/%d\n\n',kt,mtsim)
    
    F=fd(kt,:);
    R=fo(kt,:);

    % identify optimal relocation flux
    x=optimalrelocationfluxes(F,R,T,60,approx);

    % read results
    [Fs,Rs,Vr]=find(x);

    % distance of relocation
    dkemd(kt)=sum(T(sub2ind(size(T),Fs,Rs)).*Vr);

end

end

