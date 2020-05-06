


function Sim=calcwaitingtimes(P,c,dfinal)

n=sqrt(size(dfinal,1));
tsim=size(dfinal,2)-1;

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
dfinal=[dfinal zeros(n^2,1)];

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

            davantiinfila=dfinal(k,:)-(cumsum(cfinal(k,:))-num+1);
            waitingtimesteps(totalenodi(k)+num,1)=arrivi(k2); % time step of request
            waitingtimesteps(totalenodi(k)+num,2)=k;         % origin/destination pair
            waitingtimesteps(totalenodi(k)+num,3)=find((davantiinfila(arrivi(k2):end)<0),1)-1; % waiting time

            num=num+1;
            
        end
    end
end

% waiting time for each passenger in minutes
waiting=[waitingtimesteps(:,1:2) , waitingtimesteps(:,3)*P.e];

% moving average [what???]
binsize=round(10/P.e);
halfbin=floor(binsize/2);
waitingprof4=zeros(tsim,1);
for k2=1:tsim
    waitingprof4(k2)=mean(    waiting(   logical((waiting(:,1)>=k2-halfbin).*(waiting(:,1)<=k2+halfbin))   ,3)   ,'omitnan');
end
waitingprof4(isnan(waitingprof4))=0; % for bins without arrivals, assume 0 waiting


Sim.waiting=waiting;
Sim.waitingMAV10min=waitingprof4;
Sim.waitsummary=[  prctile(waiting(:,3),[0,2.5,25,50,75,97.5,100])' ; mean(waiting(:,3)) ; std(waiting(:,3))  ];
Sim.waitsummaryLgnd=["min";"2.5 pctile";"25 pctile";"50 pctile";"75 pctile";"97.5 pctile";"max";"mean";"st.dev."];


return


dfinal=X(1:n^2,:);