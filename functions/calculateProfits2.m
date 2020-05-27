function [R,Rcost,U]=calculateProfits2(Pmat,Resmat,subsidy,penalty)

DroppedUtility=-15;

R=zeros(size(Resmat));
Rcost=zeros(size(Resmat));
U=zeros(size(Resmat));

for d1=1:size(Resmat,1)
    for d2=1:size(Resmat,2)
        pooledtrips=Resmat(d1,d2).Trips.pooledtrips*Resmat(d1,d2).Trips.totalrequests;
        droppedtrips=Resmat(d1,d2).dropped*Resmat(d1,d2).Trips.totalrequests;
        nonpooled=Resmat(d1,d2).Trips.totalrequests-pooledtrips-droppedtrips;
        [~,~,R(d1,d2)]=calculateCosts3(Pmat(d1,d2),Resmat(d1,d2),0.1);
        R(d1,d2)=R(d1,d2)+(nonpooled+pooledtrips)*subsidy-droppedtrips*penalty;
        Rcost(d1,d2)=(nonpooled+pooledtrips)*subsidy-droppedtrips*penalty;
        
        U(d1,d2)=sum(Resmat(d1,d2).Economics.utility(logical(1-full(Resmat(d1,d2).Trips.dropped))))+sum(Resmat(d1,d2).Trips.dropped)*DroppedUtility;
        
    end
end