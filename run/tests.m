%% check for constraint violations

d=Res.Internals.d;

% charging when not in charging stations and with delay=0
c(1)=sum(sum((d>0).*Res.Internals.s2));

% charging more than possible
c(2)=sum(abs(round(Res.Sim.e(:),5))>P.Tech.chargekw);

% charging with wrong status
c(3)=sum(sum((abs(round(Res.Sim.e,5))>0).*(Res.Internals.s2(1:end-1,:)==0)));

% double status
c(4)=sum(sum((Res.Internals.s1+Res.Internals.s2+Res.Internals.s3)>1));

% moving while charging
c(5)=sum(sum((d>0).*Res.Internals.s2));

% soc limits
c(6)=(max(Res.Sim.q(:))>P.Operations.maxsoc)+(min(Res.Sim.q(:))<P.Operations.minsoc);

% 
% c(7)=Res

% delay change should be not more than 1
c(8)=sum((max(Res.Internals.d(1:end-1,:)-Res.Internals.d(2:end,:)))>1);

if sum(c)>0
    c
    warning('constraint violations!');
else
    fprintf('no violations\n');
end

if 0

%% specific trip

v=14;
z=900:1000;
% d s1 s2 s3
full([Res.Internals.d(z,v) Res.Internals.s1(z,v) Res.Internals.s2(z,v) Res.Internals.s3(z,v)])
plot(Res.Internals.d(z,v))

end




