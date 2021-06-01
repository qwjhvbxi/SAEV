function [fo,fd]=computeod(A,Atimes,cumulativeTripArrivals,n,Beta)

if nargin<3 || isempty(cumulativeTripArrivals)
    % fast method for trip generation
    maxt=max(Atimes(:));
    Abuck=histc(Atimes(:,1),1:maxt+1); % number of trips in each minute
    cumulativeTripArrivals=[0;cumsum(Abuck)]; % total number of trips up to each minute
end
if nargin<4 || isempty(n)
    n=max(max(A(:,1:2)));
end
if nargin<5 || isempty(Beta)
    Beta=1;
end

t=max(Atimes(:,1));
fo=zeros(t,n);
fd=zeros(t,n);

for kt=1:t
    
    thisStep=cumulativeTripArrivals((kt-1)*Beta+1)+1:cumulativeTripArrivals(min(length(cumulativeTripArrivals),kt*Beta+1));
    thisNumber=length(thisStep);
    
    % feeder and receiver stations
    fo(kt,:)=accumarray([A(thisStep,1);n],[ones(thisNumber,1);0])';
    fd(kt,:)=accumarray([A(thisStep,2);n],[ones(thisNumber,1);0])';
    
end