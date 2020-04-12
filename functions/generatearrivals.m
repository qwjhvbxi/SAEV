% input: P.n, P.tsim, P.tripsperhour; output: passengers

% generate arrivals
passengers=zeros(P.n,P.n,tsim+100);
factorialvec=factorial(0:60);
freq=tripsperhour/60*P.e;
N=0:max(10,max(freq)*3);
for j3=1:tsim+100

    hh=rem(ceil(j3/60*P.e)-1,24)+1;

    if freq(hh)<=20
        poissonCMF=[cumsum(freq(hh).^N./factorialvec(1:N(end)+1).*exp(-freq(hh))) , 1];
        expect=find(poissonCMF>rand(),1)-1;
    else
        expect=max(0,round(freq(hh)+randn()*sqrt(freq(hh))));
    end

    % extract indices of nodes for departure and arrival point
    [Q]=discretesample(tripsClusteredInterzonal{hh}(:,1),expect);
    nodesind=sub2ind([P.n,P.n],tripsClusteredInterzonal{hh}(Q,2),tripsClusteredInterzonal{hh}(Q,3));

    passengers(:,:,j3)=reshape(histc(nodesind,1:P.n^2),P.n,P.n);

end