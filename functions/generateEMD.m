%% Trips=generateEMD(A,Atimes,T,etsim,FileName)
% dkemd, dkod, dktrip are the number of minutes of travel for
% relocation, serving trips, and total, respectively, for each
% energy layer time step. 
%
% see also generalC

function Trips=generateEMD(A,Atimes,T,Beta,FileName,chargingStations,Clusters)

DataFolder=setDataFolder();

if nargin==7
    As=Clusters(A);
    Ts=T(chargingStations,chargingStations);
else
    As=A;
    Ts=T;
end


emdname=[DataFolder 'temp/emd-' FileName '-' num2str(Beta) '.mat'];
if exist(emdname,'file')
    load(emdname,'dkemd','dkod','dktrip');
else

    % is a probability distribution of trips available?
    probabilistic=false;

    if probabilistic
        % calculate from known distribution
        error('not implemented');
    else
        dkemd=ApproxRelocDist(As,Atimes,Ts,Beta);
    end
    
    % calculate distance
    t=double(max(Atimes(:)));
    mtsim=round(t/Beta);
    t=Beta*mtsim;
    dk=zeros(t,1);
    for i=1:t
        ThisMinute=logical(Atimes(:,1)==i);
        dk(i)=sum(T(sub2ind(size(T),A(ThisMinute,1),A(ThisMinute,2))));
    end
    dkod=sum(reshape(dk(1:t),Beta,mtsim))';
    dktrip=dkemd+dkod;
    
    save(emdname,'dkemd','dkod','dktrip');
    
end

Trips.dkemd=dkemd;
Trips.dkod=dkod;
Trips.dktrip=dktrip;

end

