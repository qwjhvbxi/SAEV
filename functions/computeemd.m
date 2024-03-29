%% dkemd=COMPUTEEMD(fo,fd,T,Beta,relocationNodes)
% Compute travel time for relocation during time horizon.
% 'fo' and 'fd' are respectively the number of departures and arrivals at
% each node at each time interval. Size: [minutes x no. nodes]
% 'T' is the distance matrix.
% 'Beta' is the length (in minutes) of one interval in the horizon.
% 'relocationNodes' is the index of the nodes which are relocation
% destinations (associated with A).
% 
% See also: mainsim, computetraveltime

function dkemd=computeemd(fo,fd,T,Beta)

if ~isempty(Beta) && Beta>0

    fprintf('\n Calculating approximate relocation distance\n\n')

    [thisT]=gettraveltimenow(T,0);
    nr=size(thisT,1);
    approx=(nr>100);

    t=size(fo,1);
    mtsim=round(t/Beta);

    dkemd=zeros(mtsim,1);

    for kt=1:mtsim

        % progress report
        clc
        fprintf('\n %d/%d\n\n',kt,mtsim)

        [thisT]=gettraveltimenow(T,kt*Beta);

        thisStep=(kt-1)*Beta+1:min(t,kt*Beta);
        F=sum(fd(thisStep,:));
        R=sum(fo(thisStep,:));

        % identify optimal relocation flux
        x=optimalrelocationfluxes(F,R,thisT,60,approx);

        % read results
        [Fs,Rs,Vr]=find(x);

        % distance of relocation
        dkemd(kt)=sum(thisT(sub2ind(size(thisT),Fs,Rs)).*Vr);

    end

else 
    
    dkemd=[];

end

end

