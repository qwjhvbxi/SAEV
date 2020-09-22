function [FailMinutes,DeltaPower]=testFCR(P,R)
DataFolder=setDataFolder();
load([DataFolder 'grid/' P.FCR.filename],'f');
ReshapeFactor=size(f,1)/1440*P.e;
f=average2(f(:,P.gridday),ReshapeFactor);
FCRe=min(1,max(-1,(1-(f-P.FCR.limits(1))/(P.FCR.limits(2)-P.FCR.limits(1))*2))); % needed FCR
FCRe=FCRe(1:1440);
PowerExpect=FCRe*P.FCR.contracted*1000;
Power=sum(R.Sim.ef,2);
DeltaPower=sparse(double(round(Power-PowerExpect,7)));
FailMinutes=(sum(abs(DeltaPower)>0));