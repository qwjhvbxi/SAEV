%% dkod=COMPUTETRAVELTIME(A,Atimes,T,Beta)
% Calculate distance to travel during horizon. 
% 'Beta' is the length (in minutes) of one interval in the horizon.
% 
% See also: mainsim, computeemd

function dkod=computetraveltime(A,Atimes,T,Beta)

% last minute
t=double(max(Atimes(:)));

% length of total interval
mtsim=round(t/Beta);
% t=Beta*mtsim;

% initialize
dk=zeros(mtsim,1);

for kt=1:mtsim
    
    % progress report
    clc
    fprintf('\n %d/%d\n\n',kt,mtsim)
    
    [thisT]=gettraveltimenow(T,kt*Beta);
    
    % find trips in this time step
    thisStep=logical((Atimes(:,1)>=(kt-1)*Beta+1).*(Atimes(:,1)<=kt*Beta));
    dk(kt)=sum(thisT(sub2ind(size(thisT),A(thisStep,1),A(thisStep,2))));
    
end

dkod=dk;

end