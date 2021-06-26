%% [V,B,tripdist,tripdistkm,relodist,relodistkm,queue]=tripassignmentsaev(Vin,Bin,Par)
% Trip assignment for SAEV simulation with active rebalancing
% 
% input
% Vin: vehicles information in the form:
% [station delay soc connected]
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
% relodist: total distance for relocation to pickup
% queue: queued passengers (IDs)
%
% See also: main

function [V,B,tripdist,tripdistkm,relodist,relodistkm,queue]=tripassignmentsaev(Vin,Bin,Par)

tripdist=0;
tripdistkm=0;
relodist=0;
relodistkm=0;
queue=zeros(100,1);
ql=0;

ad=Par.consumption/Par.battery*Par.Epsilon;    % discharge rate per time step (normalized)

% if there are trips
if ~isempty(Bin)
    
    m=size(Bin,1);
    waiting=Bin(:,3);
    chosenmode=(waiting>0);
    waitingestimated=zeros(m,1);
    modeutilities=zeros(m,2);
    dropped=zeros(m,1);
    
    % limit of assignments per minute
    v=round(size(Vin,1)/2);
    queueextended=[];
    if size(Bin,1)>v

        queueextended=(v+1:m)';
        waiting(v+1:m)=waiting(v+1:m)+Par.Epsilon;
        
        % if waiting exceeds maximum, reject directly and remove from queue
        queueextended=queueextended.*(waiting(v+1:m)<Par.maxwait);
        
        Bin=Bin(1:v,:);
        
    end
    
    %n=size(Par.Tr,1);
    
    ui=Vin(:,1);
    di=Vin(:,2);
    Used=zeros(length(ui),1);
    Connected=logical(Vin(:,4));
    
    m=size(Bin,1);
    
    distancepickup=Par.Tr(ui,Bin(:,1))';
    distancetomove=Par.Tr(sub2ind(size(Par.Tr),Bin(:,1),Bin(:,2)));
    distancetomovekm=Par.D(sub2ind(size(Par.D),Bin(:,1),Bin(:,2)));
    
    EnergyReq=(distancepickup+(distancetomove+di'))*ad+Par.minsoc;
    
    % create matrix X
    X=(distancepickup + ... distance from passenger
        ones(m,1)*( di' + ... current delay
                    0.1 + ... indicate the vehicle is available
                    Vin(:,4)'*0.25*Par.chargepenalty + ... penalty for currently charging vehicles
                    (1-Vin(:,3)')*0.25) ... penalty for low soc vehicles
                   )...
         .*(EnergyReq<ones(m,1)*Vin(:,3)'); % requirement for enough SOC
     
    X(X==0)=NaN;
     
    if sum(Connected)<=Par.LimitFCR
        X(:,Connected)=NaN;
    end
    
    for tripID=1:m
        
        %if sum(~isnan(X(tripID,:)))
        % best vehicle
        [Delay,uids]=min(X(tripID,:));
        
        % if there is a possible vehicle
        if ~isnan(Delay)
            WaitingTime=floor(Delay)*Par.Epsilon;
        else
            WaitingTime=NaN;
        end

        % avoid changing chosen mode after deciding
        if chosenmode(tripID)==0

            if Par.modechoice

                if isnan(WaitingTime)
                    
                    AcceptProbability=0;
                    
                else
                    
                    Tariff=Bin(tripID,4);
                    UtilitySAEV=-Tariff-(WaitingTime+distancetomove(tripID))*Par.VOT/60;
                    AcceptProbability=exp(UtilitySAEV)/(exp(UtilitySAEV)+Bin(tripID,5));
                    
                end
                
                modeutilities(tripID,:)=[UtilitySAEV , log(Bin(tripID,5))];
                
            else

                AcceptProbability=1;

            end

            chosenmode(tripID)=(rand()<AcceptProbability);

            waitingestimated(tripID)=WaitingTime;
            
        end


        if chosenmode(tripID)==1
            
            % if the best vehicle is at the station
            if (waiting(tripID)<=Par.maxwait && isnan(WaitingTime)) || (waiting(tripID)+WaitingTime<Par.maxwait)
            % if (waiting(tripID)+WaitingTime)<=Par.maxwait
            
                if ~isnan(WaitingTime)

                    waiting(tripID)=waiting(tripID)+WaitingTime;

                    % update travelled distance
                    tripdist=tripdist+distancetomove(tripID);
                    tripdistkm=tripdistkm+distancetomovekm(tripID);

                    % update additional travel distance to pickup
                    pickupdist=Par.Tr(ui(uids),Bin(tripID,1));
                    pickupdistkm=Par.D(ui(uids),Bin(tripID,1));
                    relodist=relodist+pickupdist;
                    relodistkm=relodistkm+pickupdistkm;

                    % accept request and update vehicle position
                    ui(uids)=Bin(tripID,2); % position
                    di(uids)=di(uids)+pickupdist+distancetomove(tripID); % delay

                    % update X
                    X(:,uids)=(Par.Tr(Bin(:,1),ui(uids))+ ... 
                            di(uids)+... updated delay
                            0.1 + ... indicate the vehicle is available
                            (1-Vin(uids,3))*0.25).* ... penalty for low soc vehicles
                            (EnergyReq(:,uids)+(ad*(pickupdist+distancetomove(tripID)))<Vin(uids,3)');
                    X(X(:,uids)==0,uids)=NaN;

                    % remove vehicles that are needed for FCR
                    Connected(uids)=0;
                    if sum(Connected)<=Par.LimitFCR
                        X(:,Connected)=NaN;
                    end

                    Used(uids)=1;

                else

                    % increase waiting time (minutes) for this trip
                    waiting(tripID)=waiting(tripID)+Par.Epsilon;

                    % add this trip to the queue
                    ql=ql+1;  % current trip
                    queue(ql)=tripID;

                end
                
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
        else

            % walking

        end
        
    end
    
    queue=[queue(queue>0);queueextended(queueextended>0)];
    
    
    %% report results
    
    V=[ui , di , Used];
    B=[chosenmode , waiting , dropped , waitingestimated , modeutilities];
    
else
    
    V=Vin(:,1:2);
    B=[];
    
end


