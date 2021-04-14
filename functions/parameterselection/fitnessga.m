function x=fitnessga(tve)

load('tempP.mat','P');

P.tx=round(tve(1));
P.ts=round(tve(2));
P.tr=round(tve(3));
P.m=round(tve(4));       % number of vehicles

% launch simulation
Res=main(P,0,0);
 
x=Res.avgwait*10+Res.peakwait;


