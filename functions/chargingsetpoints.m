function [e,ef]=chargingsetpoints(Par,qi,s,zmacro,fi,newSetPoint)

if Par.fcr  % setpoint based

    if newSetPoint
        [setPoints]=setpointfleet(Par,qi,s,zmacro(1:2));
    end

    [e,ef]=setpointvehicle(Par,qi,s,fi,setPoints);


else     % capacity based
    
    m=length(qi);

    v2gallowed=qi>Par.v2gminsoc;
    extracharge=(qi<Par.refillmaxsoc);
    chargevector=max(-1,min(1,(ones(1,m)*(zmacro(1)/zmacro(3))-v2gallowed*(zmacro(2)/zmacro(3))+extracharge)))*Par.ac;

    capUp=s.*min(Par.ac,Par.maxsoc-qi); % charge
    capDown=s.*min(Par.ac,(qi-Par.minsoc)*Par.efficiency); % discharge

    e=min(capUp,max(0,chargevector))+max(-capDown,min(0,chargevector));
    ef=0;

end

end


% function [SetPoints]=setpointfleet(Par,q,s,z)
% 
% % power exchanged for vehicles charging
% acv=(q<Par.fastchargesoc)*Par.ac+(q>=Par.fastchargesoc)*Par.ac*Par.slowchargeratio;
% 
% % aggregate set point (kWh)
% SetPointUpPeriod=(z(1));  % set point of aggregate fleet (kWh)
% SetPointDownPeriod=(z(2));  % set point of aggregate fleet (kWh)
% 
% % expected total vehicle capacity in the period (kWh)
% CapUpPeriod=max(0,s.*min(acv*Par.H,Par.maxsoc-q)*Par.battery); % charge
% CapDownPeriod=max(0,s.*min(acv*Par.H,(q-Par.v2gminsoc)*Par.efficiency)*Par.battery); % discharge
% 
% % set point for each vehicle for each time step (kWh)
% SetPoints(1)=min(SetPointUpPeriod,sum(CapUpPeriod))/Par.H;  % set point of aggregate fleet (kWh) UP
% SetPoints(2)=min(SetPointDownPeriod,sum(CapDownPeriod))/Par.H;  % set point of aggregate fleet (kWh) DOWN
% 
% end
% 
% 
% function [e,ef]=setpointvehicle(Par,q,s,f,SetPoints)
% 
% % power exchanged for vehicles charging
% acv=(q<Par.fastchargesoc)*Par.ac+(q>=Par.fastchargesoc)*Par.ac*Par.slowchargeratio;
% 
% % available power from fleet
% maxsoceff=1; % maxsoceff=Par.maxsoc;
% v2gallowed=q>Par.v2gminsoc;
% 
% % actual capacity in this time step
% CapUp=s.*min(acv,maxsoceff-q); % charge
% CapDown=s.*v2gallowed.*min(acv,(q-Par.minsoc)*Par.efficiency); % discharge
% 
% % calculate ratios
% if sum(CapUp)>0
%     eRatioUp=min(1,SetPoints(1)/(sum(CapUp)*Par.battery));
% else
%     eRatioUp=0;
% end
% if sum(CapDown)>0
%     eRatioDown=min(1,SetPoints(2)/(sum(CapDown)*Par.battery));
% else
%     eRatioDown=0;
% end
% 
% % calculate charging for each vehicle
% e=CapUp*eRatioUp-CapDown*eRatioDown;
% 
% 
% %% FCR provision
% 
% % needed FCR power
% FCRNeed=(f-50)/(Par.limits(2)-Par.limits(1))*2;
% FCRNeedUp=Par.af*min(1,max(0,FCRNeed)); % charge
% FCRNeedDown=Par.af*min(1,max(0,-FCRNeed)); % discharge
% 
% % available power from fleet
% AvailableUp=s.*min(Par.ac-max(0,e),1-(q+e)); % charge
% AvailableDown=s.*min(Par.ac-max(0,-e),(q+e)-Par.minsoc); % discharge
% 
% % calculate ratios
% if sum(AvailableUp)>0
%     FCRRatioUp=min(1,FCRNeedUp/sum(AvailableUp));
% else
%     FCRRatioUp=0;
% end
% if sum(AvailableDown)>0
%     FCRRatioDown=min(1,FCRNeedDown/sum(AvailableDown));
% else
%     FCRRatioDown=0;
% end
% 
% % calculate FCR power contributions
% ef=AvailableUp.*FCRRatioUp-AvailableDown.*FCRRatioDown;
% 
% end