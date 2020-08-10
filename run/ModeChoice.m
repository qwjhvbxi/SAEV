

addpath plots
DataFolder=setDataFolder();



[P2,R2]=generateplotline3('NYC2016',[],'Operations.maxwait',20,'modechoice',[false true],'m',[8000 9000 10000]);

[~,R2d]=generateplotline3('NYC2016','dropped','Operations.maxwait',20,'modechoice',[false true],'m',[8000 9000 10000]);


for j=1:3
    Walking(j)=1-sum(R2(2,j).Sim.chosenmode)/length(R2(2,j).Sim.chosenmode);
    Served=((R2(2,j).Sim.chosenmode==1).*(R2(2,j).Sim.dropped==0));
end

[OperatingProfits,~,~]=calculateProfits2(P2,R2);


TotalReq=length(R2(1).Sim.waiting);
FleetCost=reshape([P2.m],size(OperatingProfits))*(40000/10/365);
DroppedCost=reshape([R2.dropped],size(OperatingProfits))*TotalReq*15;
ChargingCost=reshape([R2.cost],size(OperatingProfits));

mvec=[8000 9000 10000];

figure('Units','centimeters','Position',[10,7,10,7])
% plot(OperatingProfits')
ax = axes;
ax.ColorOrder = [0 0 1; 1 0 0];
hold on
a1=plot(mvec,(OperatingProfits-FleetCost)','^-.');
a2=plot(mvec,(OperatingProfits-FleetCost-DroppedCost)','x:');
a3=plot(mvec,(OperatingProfits-FleetCost-ChargingCost/100)','v--','LineWidth',1.5);
a4=plot(mvec,(OperatingProfits-FleetCost-DroppedCost-ChargingCost/100)','o-','LineWidth',1.5);
ylabel('$')
xlabel('fleet size')
xticks(mvec);
xlim([7500,10500])
ylim([-4,2]*10^5)
legend([a1(1),a2(1),a3(1),a4(1),a4(2)],{'fleet+movement';'fleet+movement+dropped';'total cash profit';'total profit (no mc)';'total profit (mc)'},'Location','SouthWest')
print([DataFolder 'figures/ModeChoice/profits'],'-dpng','-r300');

figure('Units','centimeters','Position',[10,7,10,7])
ax = axes;
ax.ColorOrder = [0 0 1; 1 0 0;1 0 0 ];
hold on
plot(mvec,R2d'*100,'x-','LineWidth',1.5)
ylabel('unserved %')
yyaxis right
plot(mvec,Walking*100,'s-','LineWidth',1.5)
ylabel('%')
xlabel('fleet size')
xticks(mvec);
xlim([7500,10500])
legend({'unserved (no mc)';'unserved (mc)';'walking'})
print([DataFolder 'figures/ModeChoice/walking'],'-dpng','-r300');



%% waiting time estimation vs real

R1=R2(2,3);
Selection=logical((R1.Sim.chosenmode==1).*(R1.Sim.dropped==0).*(R1.Sim.waitingestimated>0).*(R1.Sim.waiting>0));
prova=full([R1.Sim.waitingestimated(Selection) R1.Sim.waiting(Selection)]);
M=zeros(10,10);
for j=1:10
    for k=1:10
        M(j,k)=sum(prod(prova==[j k]*2,2));
    end
end
figure
imagesc(2:2:20,2:2:20,M)

