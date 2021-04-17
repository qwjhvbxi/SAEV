% [Res]=MAIN(P[,extsave,dispiter])
% Run SAEV simulation and relocation/charge oprimization.
% Vehicles start at beginning of time step, arrive at end of time step.
% 
% u: destination or position at beginning of time step
% q: SOC at beginning of time step
% e: charging energy exchanged during time step
% ef: FCR energy exchanged during time step
% b: imbalance from relocation module
% d: delay at beginning of time step
% status:
%     1 charging
%     2 idle
%     3 moving
%     4 moving for relocation
%     5 moving to charging station
%
% TODO: add charging station size
% 
% See also: getp

function [Res]=main(P,extsave,dispiter)

%% initializations and input check

addpath functions utilities 
DataFolder=getdatafolder();
if nargin<3
    dispiter=2;
    if nargin<2
        extsave=1;
    end
end


%% console call
% If called without variables, presents available options.

if nargin==0
    P=getp();
end


%% generate unique Hash for savefile

Hash=DataHash(P);
simname=[DataFolder 'out/' Hash '.mat'];
if extsave<2
    if extsave>=0 && exist(simname,'file')
        % if the file already exists, just return the previous result
        load(simname,'Res');
        return
    end
end

if isfield(P,'trlayeralg') && strcmp(P.trlayeralg,'opti')
    Res=generalmilp(P,extsave,dispiter);
    return
end

Res=simulation(P,extsave,dispiter);

% save results
if extsave>0
    save(simname,'Res','P');
end


%% end display

if dispiter<0
    fprintf('sim #%d ',-dispiter)
end
fprintf('successfully completed \n');