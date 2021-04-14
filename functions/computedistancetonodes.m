%% [DistanceToNode,NodeID]=computedistancetonodes(A,C,MaxNodes)
% calculate distance from coordinates in A to nodes in C.

function [DistanceToNode,NodeID]=computedistancetonodes(A,C,MaxNodes)
    
Dist=sqrt((A(:,1)*ones(1,size(C,1))-ones(size(A,1),1)*C(:,1)').^2+(A(:,2)*ones(1,size(C,1))-ones(size(A,1),1)*C(:,2)').^2);
[DistanceToNode,NodeID]=sort(Dist,2);

if nargin==3
    DistanceToNode=DistanceToNode(:,1:MaxNodes);
    NodeID=NodeID(:,1:MaxNodes);
end

end