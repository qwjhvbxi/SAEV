%% [A,Atimes]=CONVERTATIMES(passengers) 
% from matrix form of OD to list form
% 
% [passengers]=CONVERTATIMES(A,Atimes,n,tsim) 
% from list form to matrix form
% 

function [a,b]=convertAtimes(c,d,e,f)

if nargin<2
    
    % from list form to matrix form
    
    passengers=c;
    
    A=NaN;
    Atimes=NaN;
    
    a=A;
    b=Atimes;
    
    warning('not implemented');
    return
    
    
else
    
    % from matrix form of OD to list form
    
    A=c;
    Atimes=d;
    n=e; % number of nodes
    tsim=f; % number of time steps
    b=NaN;
    
    a=zeros(n,n,tsim);
    
    for i=1:tsim
        kt=(Atimes(:,1)==i);
        
        A2=sub2ind([n,n],A(kt,1),A(kt,2));
        a(:,:,i)=reshape(accumarray(A2,1,[n^2,1]),n,n);
    end
    
end

end