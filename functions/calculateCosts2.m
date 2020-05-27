%% [Revenue,Costs,Profits]=calculateCosts2(P,Res,tariff,travelcost)
% travelcost and operatorcosts are in $/min, 

function [Revenue,Costs,OperatingProfits]=calculateCosts2(P,Res,tariff,travelcost)

% retrieve trip matrix
tripFileLocation=[DataFolder 'trips/' P.tripfile '.mat'];
load(tripFileLocation,'A');
if iscell(A) && P.scenarioid>0 % file with multiple scenarios?
    A1=double(A{P.scenarioid});
else
    A1=double(A);
end

TraveledTrips=logical((Res.Sim.chosenmode==1).*(Res.Sim.dropped==0));

TraveledDistance=sum(Res.Params.Tr(sub2ind(size(Res.Params.Tr),A1(TraveledTrips,1),A1(TraveledTrips,2))));
RelocationDistance=sum(Res.Sim.relodist);

% calculate revenues
Revenue=TraveledDistance*tariff;

% calculate moving costs
Costs=(TraveledDistance+RelocationDistance)*travelcost;


OperatingProfits=Revenue-Costs;

