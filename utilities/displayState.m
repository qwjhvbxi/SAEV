function displayState(i,itot,dispiter,comptime,displayResolution)

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
    fprintf(char('-'*ones(1,displayResolution)));
    fprintf('\n')
    fprintf(char('*'*ones(1,floor(i/itot*displayResolution))));
    fprintf('\n\n')
elseif dispiter==2
    clc
    displayprogress(i,itot,comptime); %cputime-starttime); 
end
