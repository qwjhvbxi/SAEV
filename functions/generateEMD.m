%% [Trips,fo,fd]=generateEMD(A,Atimes,T,etsim,FileName)
% generate number of arrivals at each station
% generate EMD in case of aggregate energy layer
% calculate from expected arrivals
% etsim is the number of energy layer time steps in a day
% dkemd, dkod, dktrip are the number of minutes of travel for
% relocation, serving trips, and total, respectively, for each
% energy layer time step. fk
%
% see also generalC

function [Trips,fo,fd]=generateEMD(A,Atimes,T,etsim,FileName)

DataFolder=setDataFolder();
n=size(T,1);

statsname=[DataFolder 'temp/tripstats-' FileName '-N' num2str(n) '.mat'];
if exist(statsname,'file')
    load(statsname,'fo','fd','dk');
else
    [~,fo,fd,dk]=tripstats2(A,Atimes,T);
    save(statsname,'Atimes','fo','fd','dk');
end


emdname=[DataFolder 'temp/emd-' FileName '-' num2str(etsim) '.mat'];
if exist(emdname,'file')
    load(emdname,'dkemd','dkod','dktrip','fk');
else

    % is a probability distribution of trips available?
    probabilistic=false;

    if probabilistic
        % calculate from known distribution
        error('not implemented');
    else
        [dkemd,dkod,dktrip,fk]=generatetripdataAlt(fo,fd,dk,T,etsim);
    end
    save(emdname,'dkemd','dkod','dktrip','fk');

end

Trips.dkemd=dkemd;
Trips.dkod=dkod;
Trips.dktrip=dktrip;
Trips.fk=fk;

end

