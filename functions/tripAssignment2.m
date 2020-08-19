function [V,B,tripdist,queue]=tripAssignment2(Vin,Bin,Par)

tripdist=0;
queue=zeros(100,1);

% if there are trips
if ~isempty(Bin)
    
    WaitingCostToggle=1;
    
    n=size(Par.Tr,1);
    m=size(Bin,1);
    
    ui=Vin(:,1);
    di=Vin(:,2);
    Used=zeros(length(ui),1);
    
    waiting=Bin(:,3);
    chosenmode=(waiting>0);
    waitingestimated=zeros(m,1);
    dropped=zeros(m,1);
    
    distancetomove=Par.Tr(sub2ind(size(Par.Tr),Bin(:,1),Bin(:,2)));
    
    EnergyReq=distancetomove*Par.ad+Par.minsoc;
    
    % create matrix X
    X=(Par.Tr(Bin(:,1),ui) + ... distance from passenger
        ones(m,1)*( di' + ... current delay
                    Vin(:,4)'*0.5 + ... penalty for currently charging vehicles
                    (1-Vin(:,3)')*0.5) ... penalty for low soc vehicles
                   )...
         .*(EnergyReq<Vin(:,3)'); % requirement for enough SOC
     
     X(X==0)=NaN;
    
    for tripID=1:m
        
        % best vehicle
        [Delay,uids]=min(X(tripID,:));
        
        WaitingTime=floor(Delay)*Par.e;

        % avoid changing chosen mode after deciding
        if chosenmode(tripID)==0

            if Par.modechoice

                UtilitySAEV=-distancetomovesorted(ka)*Par.e*Bin(tripID,4)-WaitingTime*VOT/60*WaitingCostToggle;

                AcceptProbability=exp(UtilitySAEV)/(exp(UtilitySAEV)+Bin(tripID,5));

            else

                AcceptProbability=1;

            end

            chosenmode(tripID)=(rand()<AcceptProbability);

            waitingestimated(tripID)=WaitingTime;

        end


        if chosenmode(tripID)==1
            
            % if the best vehicle is at the station
            if WaitingTime<=Par.maxwait

                waiting(tripID)=WaitingTime;
                
                % accept request and update vehicle position
                ui(uids)=Bin(tripID,2); % position
                di(uids)=floor(Delay)+distancetomove(tripID); % delay

                % update travelled distance
                tripdist=tripdist+distancetomove(tripID);
                
                % update X
                X(:,uids)=Par.Tr(Bin(:,1),ui(uids))+di(uids);
                
                Used(uids)=1;

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
    
    %% alternative assignment
    
    % assignmentoptimal(X);
    
    
    %% report results
    
    V=[ui , di , Used];
    B=[chosenmode , waiting , dropped , waitingestimated ];
    
else
    
    V=Vin(:,1:2);
    B=[];
    
end