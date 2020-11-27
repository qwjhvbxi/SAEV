%% Trips=generateEMD(A,Atimes,T,etsim,FileName)
% dkemd, dkod, dktrip are the number of minutes of travel for
% relocation, serving trips, and total, respectively, for each
% energy layer time step. 
%
% see also generalC

function Trips=generateEMD(A,Atimes,T,Beta,FileName)

DataFolder=setDataFolder();

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
        [dkemd,dkod,dktrip]=generatetripdataAlt(A,Atimes,T,Beta);
    end
    save(emdname,'dkemd','dkod','dktrip');

end

Trips.dkemd=dkemd;
Trips.dkod=dkod;
Trips.dktrip=dktrip;

end

