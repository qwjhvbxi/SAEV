% A: OD matrix
% prices: price matrix

function [newA,SelectedTrips]=tripAcceptance(A,prices)

if sum(prices(:))==0
    
    newA=A;
    SelectedTrips=[];
    return
    
else
    
    if numel(prices)>1
        SelectedTrips=rand(size(A,1),1)>prices(sub2ind(size(prices),A(:,1),A(:,2)));
    else
        SelectedTrips=rand(size(A,1),1)>prices;
    end
    newA=A(SelectedTrips,:);

end

end