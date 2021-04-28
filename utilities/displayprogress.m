%% DISPLAYPROGRESS(i,itot,dispiter,comptime[,displayResolution])
% Prints progress of simulation

function displayprogress(i,itot,dispiter,comptime,displayResolution)

if nargin<4
    displayResolution=40;
end

% update state of simulation
if dispiter<0
    if i~=itot
        fprintf('No. %d: %0.1f%% \n',-dispiter,i/itot*100);
    else
        fprintf('No. %d: finished',-dispiter);
    end
elseif dispiter==1
    clc
    progressbar(i,tot,displayResolution);
    fprintf('\n\n')
elseif dispiter==2
    clc
    displayeta(i,itot,comptime); %cputime-starttime); 
end
