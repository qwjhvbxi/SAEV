% calculate distance

function dkod=computetraveltime(A,Atimes,T,Beta)

% last minute
t=double(max(Atimes(:)));

% length of total interval
mtsim=round(t/Beta);
t=Beta*mtsim;

% initialize
dk=zeros(t,1);

for i=1:t
    
    ThisMinute=logical(Atimes(:,1)==i);
    dk(i)=sum(T(sub2ind(size(T),A(ThisMinute,1),A(ThisMinute,2))));
end

dkod=sum(reshape(dk(1:t),Beta,mtsim))';

end