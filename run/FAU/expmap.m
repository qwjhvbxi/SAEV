

vehiclecost=15:5:40;
fcrcomp=0:50:400;

VariableContract

Cost1=sum(S(4).cost)+4500*vehiclecost*7+0*fcrcomp';
% Cost1=sum(S(4).cost)+4750*vehiclecost*7+0*fcrcomp';
% Cost3=sum(S(3).cost)+5000*vehiclecost*7-2*7*24*5*fcrcomp'-sum(sum(S(3).FCRenergy))/1000*50;
% Cost3=sum(S(3).cost)+5000*vehiclecost*7-2*7*24*5*fcrcomp';
Cost3=sum(S2(1).cost)+5000*vehiclecost*7-sum(VariableContract(1,:))*fcrcomp';


Format='-depsc2';
Resolution=[];

% figure('Units','centimeters','Position',[10,7,10,7])
% imagesc(vehiclecost,fcrcomp,Cost3-Cost1)
% caxis([-1.5e5,1.5e5])
% xlabel('vehicle cost (€/day)');
% ylabel('FCR hourly compensation (€/MW)')
% box on
% colorbar
% print([DataFolder 'figures/FAU/tradeoff2'],Format,Resolution);



figure('Units','centimeters','Position',[10,7,5,5])
imagesc(vehiclecost,fcrcomp,Cost3-Cost1)
contour(vehiclecost,fcrcomp,Cost3-Cost1,'ShowText','on')
caxis([-1.5e5,1.5e5])
xlabel('vehicle cost (€/day)');
ylabel('FCR daily price (€/MW)')
box on
print([DataFolder 'figures/FAU/tradeoff_4500'],Format,Resolution);
% print([DataFolder 'figures/FAU/tradeoff_4750'],Format,Resolution);