function [fo,fd]=computeod(A,Atimes,n)

if nargin<3
    n=max(max(A(:,1:2)));
end

t=max(Atimes(:,1));
fo=zeros(t,n);
fd=zeros(t,n);
Beta=1;

for kt=1:t
    
    % find trips in this time step
    thisStep=logical((Atimes(:,1)>=(kt-1)*Beta+1).*(Atimes(:,1)<=kt*Beta));
    
    % feeder and receiver stations
    fo(kt,:)=accumarray([A(thisStep,1);n],[ones(sum(thisStep),1);0])';
    fd(kt,:)=accumarray([A(thisStep,2);n],[ones(sum(thisStep),1);0])';
    
end