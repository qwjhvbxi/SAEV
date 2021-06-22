%% SAVEP(P[,pname])
% Save input parameter struct in folder 'data/par/'
% 
% See also: main, getp

function savep(P,pname)

addpath functions utilities
DataFolder=getdatafolder();

if nargin<2
    pname=P.scenario;
end

Paddr=[DataFolder 'par/' pname];

save(Paddr,'P');

end