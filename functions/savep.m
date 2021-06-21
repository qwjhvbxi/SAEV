function savep(P,pname)

addpath functions utilities
DataFolder=getdatafolder();

if nargin<2
    pname=P.scenario;
end

Paddr=[DataFolder 'par/' pname];

save(Paddr,'P');

% if exist(Paddr,'file')
%     load(Paddr,'P');
% end

end