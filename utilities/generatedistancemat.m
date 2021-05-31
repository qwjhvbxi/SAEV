%% Generate distance matrix from coordinates
% D=generatedistancemat(C) 

function D=generatedistancemat(C)

D=sqrt((C(:,1)-C(:,1)').^2+(C(:,2)-C(:,2)').^2);

end

% % test
% a=[0,0;1,1;-1,1];
% generatedistancemat(a)