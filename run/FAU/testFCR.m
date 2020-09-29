function [FailMinutes,DeltaPower]=testFCR(P,R)

FileName=P.FCR.filename;
GridDay=P.gridday;
Limits=P.FCR.limits;
Contracted=P.FCR.contracted;
Power=sum(R.Sim.ef,2);

DataFolder=setDataFolder();
load([DataFolder 'grid/' FileName],'f');
ReshapeFactor=size(f,1)/(length(Power));
F=average2(f(:,GridDay),ReshapeFactor);

FCRe=min(1,max(-1,(F-50)/(Limits(2)-Limits(1))*2)); % needed FCR (EV POV)
FCRe=FCRe(1:1440);
PowerExpect=FCRe*Contracted*1000;
DeltaPower=sparse(double(round(Power-PowerExpect)));
FailMinutes=(sum(abs(DeltaPower)>0));

if 1
figure
hold on 
plot(PowerExpect)
plot(Power)
end