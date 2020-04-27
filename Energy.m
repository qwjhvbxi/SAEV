

cpar
Res=generalC(P,-1,1);

% 
% [A,Atimes,ASortInd,AbuckC, ...
%     ODistToNode,ONodeID,DDistToNode,DNodeID]=...
%     generateGPStrips(P);
% 
% 
% n=size(P.T,1);
% tsim=1440;
% c=convertAtimes(A,Atimes,n,tsim);
% 
% C=reshape(c(1:n^2*tsim),n^2,tsim); % arrivals actual
% 
% 
% tic
% [Atimes,fo1,fd1]=tripstats(A,Atimes,0);
% toc
% 
% tic     
% [Atimes,fo2,fd2,dk]=tripstats2(A,Atimes,P.T);
% toc
