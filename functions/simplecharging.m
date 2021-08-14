
function e=simplecharging(Par,q,s,Z)

if ~isfield(Par,'refillmaxsoc')
    Par.refillmaxsoc=0;
end

Par.ac=Par.chargekw/Par.battery/60*Par.Epsilon;    % charge rate per time step (normalized)
% Par.ad=Par.consumption/Par.battery*Par.Epsilon;    % discharge rate per time step (normalized)

m=length(q);

v2gallowed=q>Par.v2gminsoc;
extracharge=(q<Par.refillmaxsoc);
chargevector=max(-1,min(1,(ones(1,m)*(Z(1)/Z(3))-v2gallowed*(Z(2)/Z(3))+extracharge)))*Par.ac;

capUp=s.*max(0,min(Par.ac,Par.maxsoc-q)); % charge
capDown=s.*max(0,min(Par.ac,(q-Par.minsoc)*Par.efficiency)); % discharge

e=min(capUp,max(0,chargevector))+max(-capDown,min(0,chargevector));

end