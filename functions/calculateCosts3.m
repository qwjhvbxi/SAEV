%% [Revenue,Costs,Profits]=calculateCosts2(P,Res,tariff,travelcost)
% travelcost and operatorcosts are in $/min, 

function [Revenue,Costs,OperatingProfits]=calculateCosts3(P,Res,tariff,travelcost)

DataFolder=setDataFolder();

[A,~,ASortInd,~,~,~,~,~,~]=generateGPStrips(P);



TraveledTrips=logical(full((Res.Sim.chosenmode==1).*(Res.Sim.dropped==0)));

TraveledDistance=sum(Res.Params.Tr(sub2ind(size(Res.Params.Tr),A1(TraveledTrips,1),A1(TraveledTrips,2))));
RelocationDistance=sum(Res.Sim.relodist);

% calculate revenues
Revenue=TraveledDistance*tariff;

% calculate moving costs
Costs=(TraveledDistance+RelocationDistance)*travelcost;


OperatingProfits=Revenue-Costs;

