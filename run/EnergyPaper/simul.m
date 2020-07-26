%% simulations over long period
% sims in two scenarios (NY and Germany)
% 1. without optimized charging
% 2. with optimized charging
% 3. with optimized charging and carbon pricing
% (6x simulations)

% use carbon emissions from NYISO 2018, Germany 2019 (more renewables),
% carbon price... ?

if 0
    %% prova
    Period=4:8;
    P=cpar('NYC2016');
    P.tripfolder='NYC2016'; % trips in 2019
    P.m=13000;
    P.gridfile='NY_DA_2016';
    gridoffset=0; % same year
    P.carbonprice=0;
    P.Operations.maxwait=Inf;
    % with optimized charging
    P.enlayeralg='aggregate';
    [S0a,R0a]=multiDaySim(Period,P,gridoffset);
end

%% initializations

% period
% Period=4:31;
Period=18:31+14;
gridfile_a='NY_DA_2016';
gridoffset_a=0; % same year
gridfile_b='Germany_DA_2019';
gridoffset_b=-11;%3; % 2016/01/01: Fri; 2019/01/01: Tue; 
gridfile_c='Germany_DA_2019';
gridoffset_c=-11+7*20;%3+7*20; % 2016/01/01: Fri; 2019/01/01: Tue; 

% parameters
P=cpar('NYC2016');
P.Operations.maxwait=Inf; %20
P.tripfolder='NYC2016'; % trips in 2019
P.m=13000;


%% new york grid

P.gridfile=gridfile_a;
P.carbonprice=0;

% no optimization
P.enlayeralg='no';
[S1a,R1a]=multiDaySim(Period,P,gridoffset_a);

% with optimized charging
P.enlayeralg='aggregate';
[S2a,R2a]=multiDaySim(Period,P,gridoffset_a);

% with optimized charging and carbon pricing
P.carbonprice=50; % [$/ton]
[S3a,R3a]=multiDaySim(Period,P,gridoffset_a);



%% germany grid

P.gridfile=gridfile_b;
P.carbonprice=0;

% no optimization 
P.enlayeralg='no';
[S1b,R1b]=multiDaySim(Period,P,gridoffset_b);

% with optimized charging
P.enlayeralg='aggregate';
[S2b,R2b]=multiDaySim(Period,P,gridoffset_b);

% with optimized charging and carbon pricing
P.carbonprice=50; % [$/ton]
[S3b,R3b]=multiDaySim(Period,P,gridoffset_b);




%% germany grid (summer)

P.gridfile=gridfile_c;
P.carbonprice=0;

% no optimization 
P.enlayeralg='no';
[S1c,R1c]=multiDaySim(Period,P,gridoffset_c);

% with optimized charging
P.enlayeralg='aggregate';
[S2c,R2c]=multiDaySim(Period,P,gridoffset_c);

% with optimized charging and carbon pricing
P.carbonprice=50; % [$/ton]
[S3c,R3c]=multiDaySim(Period,P,gridoffset_c);



return

%% plots

DataFolder=setDataFolder();

load(['data/eleprices/' gridfile_a '.mat'],'x','y');
Ele_a=x;
Emi_a=y;
load(['data/eleprices/' gridfile_b '.mat'],'x','y');
Ele_b=x;
Emi_b=y;
load(['data/eleprices/' gridfile_c '.mat'],'x','y');
Ele_c=x;
Emi_c=y;

figure('Units','centimeters','Position',[10,7,10,7])
hold on
box on
x=linspace(0,4*7,length(Period)*48);
plot(x,reshape(Ele_a(:,Period+gridoffset_a),length(Period)*48,1))
plot(x,reshape(Ele_b(:,Period+gridoffset_b),length(Period)*48,1))
plot(x,reshape(Ele_b(:,Period+gridoffset_c),length(Period)*48,1))
xlim([0,4*7]);
xlabel('day')
ylabel('price ($/MWh)')
legend({'NYISO','DE-W','DE-S'},'Orientation','horizontal')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/Energy/prices.eps'],'-depsc2');

figure('Units','centimeters','Position',[10,7,10,7])
hold on
box on
x=linspace(0,4*7,length(Period)*48);
plot(x,reshape(Emi_a(:,Period+gridoffset_a),length(Period)*48,1))
plot(x,reshape(Emi_b(:,Period+gridoffset_b),length(Period)*48,1))
plot(x,reshape(Emi_b(:,Period+gridoffset_c),length(Period)*48,1))
xlim([0,4*7]);
xlabel('day')
ylabel('carbon intensity (g/kWh)')
legend({'NYISO','DE-W','DE-S'},'Orientation','horizontal')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/Energy/carbonintensity.eps'],'-depsc2');

% figure('Units','centimeters','Position',[10,7,10,7])
% hold on
% plot(S1.cost+S1.emissions*50,'x-')
% plot(S2.cost+S2.emissions*50,'s-')
% plot(S3.cost,'o-')


% figure('Units','centimeters','Position',[10,7,10,7])
% hold on
% x=1:3;
% plot(x,[sum(S1.cost) , sum(S2.cost), sum(S3.cost-S3.emissions*50)]/length(Period),'s:','LineWidth',1.5)
% plot(x,[sum(S1.cost+S1.emissions*50) , sum(S2.cost+S2.emissions*50), sum(S3.cost)]/length(Period),'x:','LineWidth',1.5)
% yyaxis right
% plot(x,[sum(S1.emissions) , sum(S2.emissions), sum(S3.emissions)]/length(Period),'o:','LineWidth',1.5)
% xlim([0.5,3.5]);
% xticks(x);
% xticklabels({'B','O','O+CP'})
% legend({''})



% 
% sum(S2.emissions)
% sum(S3.emissions)
% 
% figure('Units','centimeters','Position',[10,7,10,7])
% hold on
% plot(S1.avgwait)
% plot(S2.avgwait)
% plot(S3.avgwait)
% 
% figure('Units','centimeters','Position',[10,7,10,7])
% hold on
% plot(S1.dropped)
% plot(S2.dropped)
% plot(S3.dropped)
% 
% figure
% hold on
% plot(S1.emissions)
% plot(S2.emissions)
% plot(S3.emissions)
% 
% 
% 
% figure
% hold on
% plot(sum(R2(2).Sim.e,2))
% plot(sum(R3(2).Sim.e,2))
% yyaxis right
% plot(R3(2).Params.co2)



%% table

writeLine=@(S,Tax) fprintf('& %d & %d & %.0f & %.3f & %.0f & %.2f \\\\\n',full(S.totwaitprctile(3)) , full(S.totwaitprctile(4)) , mean(S.cost+S.emissions*Tax) , sum(S.cost+S.emissions*Tax)/S.totminutes*60 , mean(S.emissions) , sum(S.emissions)/S.totminutes*10^6);

writeLine(S0a,0)

% NY no carbon price
writeLine(S1a,0)
writeLine(S2a,0)

% germany winter no carbon price
writeLine(S1b,0)
writeLine(S2b,0)

% germany summer no carbon price
writeLine(S1c,0)
writeLine(S2c,0)


% NY $50 carbon tax
writeLine(S1a,50)
writeLine(S3a,0)

% germany winter $50 carbon tax
writeLine(S1b,50)
writeLine(S3b,0)

% germany summer $50 carbon tax
writeLine(S1c,50)
writeLine(S3c,0)




return


%% plot for paper

params
P.trlayeralg='simplified';
P.tx=1;
P.ts=2;
P.tr=5;
P.bmin=0;
P.m=10000;
P.scenario='NY';
load(['data/scenarios/' P.scenario '.mat'],'T');
P.T=round(T);
Res4=general11(P,-1,0);


%%

tsim=720;
x=linspace(0,24,tsim);
x1=linspace(0,24,tsim+1);
x2=linspace(0,24,1440);
xt=0:4:24;

% power exchanged
figure('Units','centimeters','Position',[10,7,10,7])
plot(x,sum(Res4.Sim.e,2)/1000)
ylabel('power (MW)')
hold on
yyaxis right
plot(x,Res4.Params.elep(1:tsim),'r-')
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('electricity price (yen/kWh)')
legend({'power','price'},'Orientation','horizontal')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print -depsc2 figures/power.eps

% waiting 
figure('Units','centimeters','Position',[10,7,10,7])
plot(x,Res4.Sim.waitingMAV10min)
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('average wait time (min.)')
% legend({'power','price'},'Orientation','horizontal','Location','best')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
print -depsc2 figures/waiting.eps

% soc 
figure('Units','centimeters','Position',[10,7,10,7])
plot(x1,mean(Res4.Sim.q,2))
xlim([0,24])
xticks(xt)
xlabel('hours')
ylabel('average SOC')
% legend({'power','price'},'Orientation','horizontal','Location','best')
set(gca,'FontUnits','points','FontWeight','normal','FontName','Times');
print -depsc2 figures/soc.eps


%% load trips only

[A,Atimes,ASortInd,AbuckC,Distances]=generateGPStrips(P);



