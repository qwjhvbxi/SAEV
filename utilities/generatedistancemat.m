%% Generate distance matrix from coordinates
% D=generatedistancemat(C) 
% D=generatedistancemat(C1,C2)
% C is a n x 2 matrix with Euclidean coordinates

function D=generatedistancemat(C1,C2)

if nargin<2
    C2=C1;
end

D=sqrt((C1(:,1)-C2(:,1)').^2+(C1(:,2)-C2(:,2)').^2);

end

% % test
% a=[0,0;1,1;-1,1];
% generatedistancemat(a)