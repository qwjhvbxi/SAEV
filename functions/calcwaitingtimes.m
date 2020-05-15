%% [waiting,tripinfo]=calcwaitingtimes(e,c,d)
% calculate waiting time for each passenger from aggregate matrices for
% optimal transport layer algorithm. 'e' is the step length in minutes; 'c'
% is the passenger arrivals, 'd' is the passengers waiting at stations

function [waiting,tripinfo]=calcwaitingtimes(e,c,d)

n=sqrt(size(d,1));
tsim=size(d,2)-1;

if min(size(c))==1
    cfinal=reshape(c(1:n^2*tsim),n^2,tsim);
else 
    if ndims(c)==3
        cfinal=reshape(c(:,:,1:tsim),n^2,tsim);
    end
end


%% calculate waiting times

% arrivals [origin/destinationPair  timesteps]
cfinal=[zeros(n^2,1) cfinal zeros(n^2,1)];

% waiting at stations [origin/destinationPair  timesteps]
d=[d zeros(n^2,1)];

% total arrivals
totalarrivals=sum(sum(cfinal));

% cumulative arrivals at each origin/destinationPair
totalenodi=[0 ; cumsum(sum(cfinal,2))];

% total waiting minutes
waitingtimesteps=zeros(totalarrivals,3);

% for each origin/destinationPair
for k=1:n^2
    
    % find moments with arrivals at origin/destinationPair
    arrivi=find(cfinal(k,:));
    
    % how many arrivals for each moment
    numeroarrivi=cfinal(k,(cfinal(k,:)>0));
    
    % initialize
    num=1;
    
    % for each moment when there are arrivals
    for k2=1:length(arrivi)

        for j=1:numeroarrivi(k2)

            davantiinfila=d(k,:)-(cumsum(cfinal(k,:))-num+1);
            waitingtimesteps(totalenodi(k)+num,1)=arrivi(k2); % time step of request
            waitingtimesteps(totalenodi(k)+num,2)=k;         % origin/destination pair
            waitingtimesteps(totalenodi(k)+num,3)=find((davantiinfila(arrivi(k2):end)<0),1)-1; % waiting time

            num=num+1;
            
        end
    end
end

% waiting time for each passenger in minutes
waiting=waitingtimesteps(:,3)*e;

% time step of request, origin/destination pair
tripinfo=waitingtimesteps(:,1:2);




