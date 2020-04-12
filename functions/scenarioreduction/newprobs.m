
% calculate new probabilities

function [q,Id,Er]=newprobs(C,p,u,nr)

n=size(C,1);
Id=zeros(n,1);
q=zeros(length(u),1);
q(nr+1:end)=NaN;

for j=1:n
    
    [a,b]=min(C(j,u(1:nr)));
    
    Id(j)=b;
    Er(j)=a;
    
    q(b)=q(b)+p(j);
    
end

% % plot similar days
% k=1
% plot(pday(:,Id==k))

end