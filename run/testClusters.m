
%% sensitivity analysis: number of clusters

% generateClusters('NYC2016',50)

Pmat{1}=cpar('NYC2016');

Pmat{2}=cpar('NYC2016-50clusters');
Pmat{3}=cpar('NYC2016-20clusters');
Pmat{4}=cpar('NYC2016-10clusters');

for i=1:length(Pmat)
    Res{i}=generalC(Pmat{i},1,2);
end

R=cell2mat(Res);
S=[R.Sim]

figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot([[R.avgwait]'],'o-')
ylabel('avg. waiting time (min)')
yyaxis right
plot([R.cost],'x-')
xticks(1:length(Pmat))
xticklabels({'200','50','20','10'})
xlabel('number of charging stations')
ylabel('charging cost')
legend({'waiting','cost'})

figure('Units','centimeters','Position',[10,7,10,7])
plot(sum([S.relodist]),'o-')
xticks(1:length(Pmat))
xticklabels({'200','50','20','10'})
xlabel('number of charging stations')
ylabel('relocation travel time (min)')
% legend({'waiting','cost'})


%% clusters + pricing

P=cpar('NYC2016-20clusters');
P.modechoice=1;

Pricing.relocationcost=0.1;
Pricing.basetariff=0.25;
Pricing.alternative=0.25;
Pricing.tp=30;
Pricing.pricingwaiting=false;
Pricing.dynamic=true;
Pricing.VOT=15;
P.Pricing=Pricing;

P.Pricing.dynamic=false;
Res1=generalC(P,-1,2);

P.Pricing.dynamic=true;
Res2=generalC(P,-1,2);



