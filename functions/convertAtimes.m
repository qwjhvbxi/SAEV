%% [A,Atimes]=CONVERTATIMES(passengers) 
% from matrix form of OD to list form
% 
% [passengers]=CONVERTATIMES(A,Atimes,n,tsim) 
% from list form to matrix form
% 

function [a,b]=convertAtimes(c,d,e)

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
    tsim=1440; % number of time steps
    b=NaN;
    
    a=zeros(n,n,tsim,'uint16');

%     for i=1:length(A)
%         a(A(i,1),A(i,2),Atimes(i,1))=a(A(i,1),A(i,2),Atimes(i,1))+1;
%     end
    
    for i=1:tsim
        kt=(Atimes(:,1)==i);
        
        A2=sub2ind([n,n],A(kt,1),A(kt,2));
        a(:,:,i)=reshape(accumarray(A2,1,[n^2,1]),n,n);
    end
    
end

end