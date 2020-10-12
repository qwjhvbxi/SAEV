%% [V,B,tripdist,relodist,queue]=tripAssignment2(Vin,Bin,Par)
% Trip assignment for SAEV simulation with active rebalancing
% 
% input
% Vin: vehicles information in the form:
% [station delay soc connected]
% Bin: passengers info in the form:
% [O D waiting]
% Par: parameters:
%     Tr, ad, e, minsoc, maxwait, modechoice 
% 
% output 
% V: vehicle movements in the form [station delay used]
% B: passengers status in the form:
% [chosenmode waiting dropped waitingestimated]
% tripdist: total distance with passengers
% relodist: total distance for relocation to pickup
% queue: queued passengers (IDs)
%
% See also: generalC

function [V,B,tripdist,relodist,queue]=tripAssignment2(Vin,Bin,Par)

tripdist=0;
relodist=0;
queue=zeros(100,1);

% if there are trips
if ~isempty(Bin)
    
    n=size(Par.Tr,1);
    m=size(Bin,1);
    v=size(Vin,1);
    
    ui=Vin(:,1);
    di=Vin(:,2);
    Used=zeros(length(ui),1);
    Connected=logical(Vin(:,4));
    if isfield(Par,'LimitFCR')
        LimitFCR=Par.LimitFCR;
    else
        LimitFCR=0;
    end
    
    waiting=Bin(:,3);
    chosenmode=(waiting>0);
    waitingestimated=zeros(m,1);
    dropped=zeros(m,1);
    
    distancepickup=Par.Tr(ui,Bin(:,1))';
    distancetomove=Par.Tr(sub2ind(size(Par.Tr),Bin(:,1),Bin(:,2)));
    
    EnergyReq=(distancepickup+(distancetomove+di'))*Par.ad+Par.minsoc;
    
    % create matrix X
    X=(distancepickup + ... distance from passenger
        ones(m,1)*( di' + ... current delay
                    Vin(:,4)'*0.5 + ... penalty for currently charging vehicles
                    (1-Vin(:,3)')*0.5) ... penalty for low soc vehicles
                   )...
         .*(EnergyReq<ones(m,1)*Vin(:,3)'); % requirement for enough SOC
     
     X(X==0)=NaN;
    
    for tripID=1:m
        
        %if sum(~isnan(X(tripID,:)))
        % best vehicle
        [Delay,uids]=min(X(tripID,:));
        
        % if there is a possible vehicle
        if ~isnan(Delay)
            WaitingTime=floor(Delay)*Par.e;
        else
            WaitingTime=NaN;
        end

        % avoid changing chosen mode after deciding
        if chosenmode(tripID)==0

            if Par.modechoice

                if isnan(WaitingTime)
                    
                    AcceptProbability=0;
                    
                else
                    
                    % Tariff=distancetomove(tripID)*Par.e*Bin(tripID,4);
                    Tariff=Bin(tripID,4);

                    UtilitySAEV=-Tariff-WaitingTime*Par.VOT/60*Par.WaitingCostToggle;

                    AcceptProbability=exp(UtilitySAEV)/(exp(UtilitySAEV)+Bin(tripID,5));
                    
                end
                
            else

                AcceptProbability=1;

            end

            chosenmode(tripID)=(rand()<AcceptProbability);

            waitingestimated(tripID)=WaitingTime;

        end


        if chosenmode(tripID)==1
            
            % if the best vehicle is at the station
            if WaitingTime<=Par.maxwait && ~isnan(WaitingTime)
                
                waiting(tripID)=WaitingTime;
                
                % update travelled distance
                tripdist=tripdist+distancetomove(tripID);
                
                % update additional travel distance to pickup
                pickupdist=Par.Tr(ui(uids),Bin(tripID,1));
                relodist=relodist+pickupdist;
                
                % accept request and update vehicle position
                ui(uids)=Bin(tripID,2); % position
                di(uids)=di(uids)+pickupdist+distancetomove(tripID); % delay
                
                % update X
                X(:,uids)=(Par.Tr(Bin(:,1),ui(uids))+di(uids)).*(EnergyReq(:,uids)+(Par.ad*(pickupdist+distancetomove(tripID)))<Vin(uids,3)');
                X(X(:,uids)==0,uids)=NaN;
                
                % remove vehicles that are needed for FCR
                Connected(uids)=0;
                if sum(Connected)<=LimitFCR
                    X(:,Connected)=NaN;
                end
                
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
    
    
    %% report results
    
    V=[ui , di , Used];
    B=[chosenmode , waiting , dropped , waitingestimated ];
    
else
    
    V=Vin(:,1:2);
    B=[];
    
end