%% [dkemd,dkod]=COMPUTEEMD(A,Atimes,T,Beta,relocationNodes)
% Compute travel time for relocation during time horizon
% 
% See also: main

function [dkemd,dkod]=computeemd(A,Atimes,T,Beta,relocationNodes)

fprintf('\n Calculating approximate relocation distance\n\n')

t=double(max(Atimes(:)));
mtsim=round(t/Beta);

dkemd=zeros(mtsim,1);
dkod=zeros(mtsim,1);

for kt=1:mtsim
    
    % progress report
    clc
    fprintf('\n %d/%d\n\n',kt,mtsim)
    
    [thisT]=gettraveltimenow(T,kt*Beta);
    
    if nargin<5
        n=size(thisT,1);
        relocationNodes=(1:n)';
    end
    thisTrelocation=thisT(relocationNodes,relocationNodes);
    nr=size(thisTrelocation,1);
    approx=(nr>100);
    
    % find trips in this time step
    thisStep=logical((Atimes(:,1)>=(kt-1)*Beta+1).*(Atimes(:,1)<=kt*Beta));
    
    % feeder and receiver stations
    F=accumarray([A(thisStep,2);nr],[ones(sum(thisStep),1);0])';
    R=accumarray([A(thisStep,1);nr],[ones(sum(thisStep),1);0])';

    % identify optimal relocation flux
    x=optimalrelocationfluxes(F,R,thisTrelocation,60,approx);

    % read results
    [Fs,Rs,Vr]=find(x);

    % distance of relocation
    dkemd(kt)=sum(thisTrelocation(sub2ind(size(thisTrelocation),Fs,Rs)).*Vr);
    
    % distance to travel
    dkod(kt)=sum(thisTrelocation(sub2ind(size(thisTrelocation),A(thisStep,1),A(thisStep,2))));

end

end

