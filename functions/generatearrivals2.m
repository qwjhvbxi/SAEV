%% [passengers]=GENERATEARRIVALS2(Trips,tph,n,tripsperday)
% Trips: cell array with [weight, origin, destination] at each hour (cell)
% tph: probability of trip starting at each hour of the day
% n: nodes
% tripsperday: number of total trips in a day

function [passengers]=generatearrivals2(OD,TPH,n,tripsperday)

tsim=1440+180;
passengers=zeros(n,n,tsim);

freqraw=interp1(1:24,TPH',linspace(1,24,1440),'pchip');
freq=repmat(freqraw/sum(freqraw)*tripsperday,1,2);
% figure
% hold on
% plot(freq)
% freq=interp1(1:24,tph',linspace(1,24,1440))*tripsperday/1440;
% plot(freq)
% freq=interp1(1:24,tph',linspace(1,24,1440),'spline')*tripsperday/1440;
% plot(freq)

N=0:max(10,max(freq)*3);
factorialvec=factorial(0:N(end)+1);%60
for j3=1:tsim
    
    % how many people in this time step?
    if freq(j3)<=20
        poissonCMF=[cumsum(freq(j3).^N./factorialvec(1:N(end)+1).*exp(-freq(j3))) , 1];
        expect=find(poissonCMF>rand(),1)-1;
    else
        expect=max(0,round(freq(j3)+randn()*sqrt(freq(j3))));
    end
    
    % where are they going?
    % extract indices of nodes for departure and arrival point
    hh=rem((rem(j3-1,60)+1)-1,24)+1;
    q=discretesample(OD{hh}(:,1),expect);
    nodesind=sub2ind([n,n],OD{hh}(q,2),OD{hh}(q,3));

    passengers(:,:,j3)=reshape(histc(nodesind,1:n^2),n,n);

end

end


