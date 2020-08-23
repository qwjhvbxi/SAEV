function tests(P,Res,v)

if nargin<=2

    DataFolder=setDataFolder();
    load([DataFolder 'scenarios/' P.scenario '.mat'],'T');
    n=size(T,1);

    if isfield(P,'clusters')
        chargingStations=P.chargingStations;
        Clusters=P.clusters;
    else
        chargingStations=(1:n)';
        Clusters=(1:n)';
    end

    %% check for constraint violations

    % charging when not in charging stations and with delay=0
    c(1)=sum(sum((Res.Internals.d>0).*Res.Internals.s2));

    % charging more than possible
    c(2)=sum(abs(round(Res.Sim.e(:),4))>P.Tech.chargekw);

    % charging with wrong status
    c(3)=sum(sum((abs(round(Res.Sim.e,5))>0).*(Res.Internals.s2(1:end-1,:)==0)));

    % double status
    c(4)=sum(sum((Res.Internals.s1+Res.Internals.s2+Res.Internals.s3)>1));

    % moving while charging
    c(5)=sum(sum((Res.Internals.d>0).*Res.Internals.s2));

    % soc limits (can go loewer if going to charging station)
    c(6)=(max(Res.Sim.q(:))>P.Operations.maxsoc)+sum(sum(Res.Sim.q<P.Operations.minsoc.*full(Res.Internals.s2+Res.Internals.s3==0)));

    % charging outside of charging stations
    c(7)=sum(setdiff(unique(full(abs(round(Res.Sim.e,5))>0).*double(Res.Sim.u(1:end-1,:))),chargingStations))>0;

    % delay change more than 1
    c(8)=sum((max(Res.Internals.d(1:end-1,:)-Res.Internals.d(2:end,:)))>1);

    % assigned distances different from traveled (not relevant with new d)
    % c(9)=(  (sum(sum(Res.Internals.d>0))+sum(max(Res.Internals.d(end,:)-1,0)))*P.e ~= (sum(Res.Sim.relodist)+sum(Res.Sim.tripdist)));

    % energy in very different from energy out
    energyIn=(full(sum(max(0,Res.Sim.e(:)))/60*P.e)+double(sum(Res.Sim.q(1,:))*P.Tech.battery));
    energyOut=(-full(sum(min(0,Res.Sim.e(:)))/60*P.e)+(sum(Res.Sim.relodist+Res.Sim.tripdist)-sum(max(Res.Internals.d(end,:),0)))*P.Tech.consumption + sum(Res.Sim.q(end,:))*P.Tech.battery);
    c(10)=(round(energyIn/energyOut,3)~=1);

    if sum(c)>0
        c
        warning('constraint violations!');
    else
        fprintf('no violations\n');
    end

else

    %% specific trip

    z=1:Res.Params.tsim;
    pos=Res.Sim.u(z,v);
    % d s1 s2 s3
    full([double(Res.Sim.q(z,v)) Res.Internals.d(z,v) Res.Internals.s1(z,v) Res.Internals.s2(z,v) Res.Internals.s3(z,v)])
    plot(Res.Internals.d(z,v))

end

end


