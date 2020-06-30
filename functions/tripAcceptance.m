% A: OD matrix
% prices: price matrix

function newA=tripAcceptance(A,prices)

if sum(prices(:))==0
    
    newA=A;
    return
    
else
    
    SelectedTrips=rand(size(A,1),1)>prices(sub2ind(size(prices),A(:,1),A(:,2)));
    newA=A(SelectedTrips,:);

end

end