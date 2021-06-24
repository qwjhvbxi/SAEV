% [Res]=MAIN(P[,extsave,dispiter])
% Run SAEV simulation and relocation/charge optimization. 
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
% See also: getp, mainsim
%
% Riccardo Iacobucci 
% https://riccardoiacobucci.com

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
    if isempty(P)
        return
    end
end


%% generate unique Hash for savefile

Warn=0;
if exist('DataHash','file')==2
    Hash=DataHash(P);
    simname=[DataFolder 'out/' Hash '.mat'];
    if extsave<2
        if extsave>=0 && exist(simname,'file')
            % if the file already exists, just return the previous result
            load(simname,'Res');
            return
        end
    end
else
    Warn=(extsave>=0);
    extsave=-1;
end


%% launch alternative benchmark simulation

if isfield(P,'trlayeralg') && strcmp(P.trlayeralg,'opti')
    Res=generalmilp(P,extsave,dispiter);
    return
end


%% default values

% TODO: define all default values here

% check optional info
if ~isfield(P,'Pricing') || isempty(P.Pricing)
    P.Pricing=struct('movingcostkm',0,'basetariffkm',0,'VOT',0,'alternativecostkm',0,'alternativecost',[],'dynamic',0,'mintariff',0);
end

if isfield(P,'FCR') && ~isempty(P.FCR) && ~isfield(P.FCR,'aggregatechargeratio')
   	P.FCR.aggregatechargeratio=1;
end


%% legacy pricing conversions

if ~isfield(P.Pricing,'mintariff') 
    P.Pricing.mintariff=0;
end

if isfield(P.Pricing,'alternative') 
    
    % convert price per minutes to prices per km;
    avgspeedkmh=30;
    mintokm=60/avgspeedkmh;
    P.Pricing.movingcostkm=P.Pricing.relocationcost*mintokm; % TODO: add km traveled!
    P.Pricing.basetariffkm=P.Pricing.basetariff*mintokm;
    
    % convert trip alternative to OD pricing
    if numel(P.Pricing.alternative)==1
        P.Pricing.alternativecostkm=P.Pricing.alternative*mintokm;
        P.Pricing.alternativecost=[];
    else
        P.Pricing.alternativecost=P.Pricing.alternative;
    end
    
    % convert dynamic/nodebased to algs;
    % P.Pricing.alg=[];
    % if P.Pricing.dynamic && ~P.Pricing.nodebased
    %     P.Pricing.alg='od';
    % end
    % if P.Pricing.dynamic && P.Pricing.nodebased
    %     P.Pricing.alg='nodes';
    % end
end


%% launch main simulation

Res=mainsim(P,dispiter);

% save results
if extsave>0
    save(simname,'Res','P');
end

% end display

if dispiter<0
    fprintf('sim #%d ',-dispiter)
end
fprintf('successfully completed \n');

if Warn
    warning('DataHash not found. Used extsave=-1.')
end

