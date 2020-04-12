function displayState(i,itot,dispiter,comptime,DisplayResolution)

if nargin<4
    DisplayResolution=40;
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
    fprintf(char('-'*ones(1,DisplayResolution)));
    fprintf('\n')
    fprintf(char('*'*ones(1,floor(i/itot*DisplayResolution))));
    fprintf('\n\n')
elseif dispiter==2
    clc
    displayprogress(i,itot,comptime); %cputime-starttime); 
end
