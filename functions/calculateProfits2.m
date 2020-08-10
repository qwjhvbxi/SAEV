%% [R,Rcost,U]=calculateProfits2(Pmat,Resmat,subsidy,penalty)
% 

function [R,Rcost,U]=calculateProfits2(Pmat,Resmat,subsidy,penalty)

DroppedUtility=-15;

R=zeros(size(Resmat));
Rcost=zeros(size(Resmat));
U=zeros(size(Resmat));

for d1=1:size(Resmat,1)
    for d2=1:size(Resmat,2)
        droppedtrips=sum(Resmat(d1,d2).Sim.dropped);
        [~,~,R(d1,d2)]=calculateCosts3(Pmat(d1,d2),Resmat(d1,d2),0.2,0.1);
        
%         Rcost(d1,d2)=(nonpooled+pooledtrips)*subsidy-droppedtrips*penalty;
%         U(d1,d2)=sum(Resmat(d1,d2).Economics.utility(logical(1-full(Resmat(d1,d2).Trips.dropped))))+sum(Resmat(d1,d2).Trips.dropped)*DroppedUtility;
        
    end
end

end



%% [Revenue,Costs,Profits]=calculateCosts2(P,Res,tariff,travelcost)
% travelcost and operatorcosts are in $/min, 

function [Revenue,Costs,OperatingProfits]=calculateCosts3(P,Res,tariff,travelcost)

[A,~,ASortInd,~,~]=generateGPStrips(P);

A=reorderVectors(A,ASortInd);

TraveledTrips=logical(full((Res.Sim.chosenmode==1).*(Res.Sim.dropped==0)));

TraveledDistance=sum(Res.Params.Tr(sub2ind(size(Res.Params.Tr),A(TraveledTrips,1),A(TraveledTrips,2))));
RelocationDistance=sum(Res.Sim.relodist);

% calculate revenues
Revenue=TraveledDistance*tariff;

% calculate moving costs
Costs=(TraveledDistance+RelocationDistance)*travelcost;


OperatingProfits=Revenue-Costs;

end