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
