%% Plot occupancy of charging stations

function csoccupancy=plotcsoccupancy(Res)

ncs=length(Res.Params.chargingStations);
csoccupancy=zeros(1440,ncs);

for i=1:1440
    mat1=logical((Res.Sim.u(i,:)==Res.Params.chargingStations(:,1)).*(Res.Internals.d(i,:)==0));
    % atChargingStation=sum(mat1);
    whichcs=(1:ncs)*mat1;
    csoccupancy(i,:)=histc(whichcs,1:ncs);
end

figure
hold on
plot(csoccupancy./Res.Params.chargingStations(:,2)')
xlim([1,1440])
line([1,1440],[1,1],'Color','k','LineStyle','--')

end

