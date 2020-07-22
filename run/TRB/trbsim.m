

P=cpar('NYC2018');
P.Operations.maxwait=20;
P.m=2000;
P.TransportLayer.tp=15;
P.enlayeralg='no';
P.TransportLayer.relocationcost=0.1;
P.TransportLayer.basetariff=0.5;
P.tripday=3;

P.pricing=false;
Res1=generalC(P,2,2)

P.pricing=true;
Res2=generalC(P,2,2)


Res1.Sim.revenues-Res1.Sim.relocationcosts
Res2.Sim.revenues-Res2.Sim.relocationcosts
(Res2.Sim.revenues-Res2.Sim.relocationcosts)/((Res1.Sim.revenues-Res1.Sim.relocationcosts))

%%

addpath plots

[P1,R1]=generateplotline3('NYC2018',[],'tripday',1:3,'Operations.maxwait',Inf,'m',5000,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',false,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.5);
[P2,R2]=generateplotline3('NYC2018',[],'tripday',1:3,'Operations.maxwait',Inf,'m',5000,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.5);

S1=[R1.Sim];
S2=[R2.Sim];

U1=sum([S1.revenues])-sum([S1.relocationcosts])
U2=sum([S2.revenues])-sum([S2.relocationcosts])

U2/U1








