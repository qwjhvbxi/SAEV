%% [V,B,tripdist,relodist,queue]=tripassignmentcarsharing(Vin,Bin,Par)
% Trip assignment for SAEV simulation
% 
% input
% Vin: vehicles information in the form:
% [station delay soc]
% Bin: passengers info in the form:
% [O D waiting]
% Par: parameters:
%     Tr, Epsilon, consumption, battery, minsoc, maxwait, modechoice 
% 
% output 
% V: vehicle movements in the form [station delay used]
% B: passengers status in the form:
% [chosenmode waiting dropped waitingestimated]
% tripdist: total distance with passengers
% relodist: total distance for relocation to pickup (0)
% queue: queued passengers (IDs)
%
% See also: main

function [V,B,tripdist,tripdistkm,relodist,queue]=tripassignmentcarsharing(Vin,Bin,Par)

tripdist=0;
tripdistkm=0;
relodist=0;
queue=zeros(100,1);

ad=Par.consumption/Par.battery*Par.Epsilon;    % discharge rate per time step (normalized)

% if there are trips
if ~isempty(Bin)
    
    WaitingCostToggle=1;
    
    % initialize current queue
    ql=0;
    n=size(Par.Tr,1);
    m=size(Bin,1);
    
    ui=Vin(:,1);
    di=Vin(:,2);
    Used=zeros(length(ui),1);
    
    waiting=Bin(:,3);
    chosenmode=(waiting>0);
    waitingestimated=zeros(m,1);
    dropped=zeros(m,1);
    
    % for each station
    for j=1:n

        % trips starting at station k
        tripsK=find(Bin(:,1)==j);

        % if there are passengers waiting at station
        if ~isempty(tripsK)

            % vehicles at this station
            uid=find(ui==j);  

            % soc of vehicles at this station (sorted by high soc and low waiting time)
            [~,usortid]=sort(Vin(uid,3)-Vin(uid,2),'descend');  % vehicle ID: uid(usortid)
            
            qj=Vin(uid(usortid),3);
            wj=Vin(uid(usortid),2); % expected waiting time to arrival

            % destination station
            destinations=Bin(tripsK,2);

            % distance of each trip
            distancetomove=Par.Tr(j,destinations);
            distancetomovekm=Par.D(j,destinations);

            % trip priority: highest waiting first, then longest travel time
            trippriority=distancetomove+waiting(tripsK)'*max(Par.Tr(:));

            % sort trips by distance (highest first)
            [~,sortid]=sort(trippriority,'descend');

            % sort trip distances
            distancetomovesorted=distancetomove(sortid);
            distancetomovekmsorted=distancetomovekm(sortid);

            % for each trip
            for ka=1:length(distancetomovesorted)

                % trip ID
                tripID=tripsK(sortid(ka));

                % candidate vehicles
                usortedi=find(qj>=distancetomovesorted(ka)*ad+Par.minsoc,1);
                uids=uid(usortid(usortedi));

                if ~isempty(uids)

                    WaitingTime=Vin(uids,2)*Par.Epsilon;

                else

                    WaitingTime=Inf;

                end


                % avoid changing chosen mode after deciding
                if chosenmode(tripID)==0

                    if Par.modechoice
                        
                        Tariff=Bin(tripID,4);
                        
                        UtilitySAEV=-Tariff-WaitingTime*Par.VOT/60*Par.WaitingCostToggle;

                        AcceptProbability=exp(UtilitySAEV)/(exp(UtilitySAEV)+Bin(tripID,5));

                    else

                        AcceptProbability=1;

                    end

                    chosenmode(tripID)=(rand()<AcceptProbability);

                    waitingestimated(tripID)=WaitingTime;

                end

                
                if chosenmode(tripID)==1

                    % if the best vehicle is at the station
                    if WaitingTime==0

                        % accept request and update vehicle position
                        ui(uids)=destinations(sortid(ka)); % position
                        di(uids)=distancetomovesorted(ka); % delay
                        qj(usortedi)=0;
                        
                        Used(uids)=1;
                        
                        % update travelled distance
                        tripdist=tripdist+distancetomovesorted(ka);
                        tripdistkm=tripdistkm+distancetomovekmsorted(ka);

                    else

                        % check if wait time of this request is less than max. wait time (minutes)
                        if waiting(tripID)<Par.maxwait

                            % increase waiting time (minutes) for this trip
                            waiting(tripID)=waiting(tripID)+Par.Epsilon;

                            % add this trip to the queue
                            ql=ql+1;  % current trip
                            queue(ql)=tripID;

                        else
                            
                            % TODO: fix pooling

                            % if max waiting exceeded, request dropped
%                             if pooling(tripID)>0
%                                 tripsdropped=find(pooling==pooling(tripID));
%                             else
                                tripsdropped=tripID;
%                             end

                            % register this as a dropped request
                            dropped(tripsdropped)=1;

                        end

                    end
                else

                    % walking

                end
            end
        end
    end
    
    V=[ui , di , Used];
    B=[chosenmode , waiting , dropped , waitingestimated ];
    
else
    
    V=Vin(:,1:2);
    B=[];
    
end