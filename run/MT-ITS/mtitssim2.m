
% P=cpar('NYC2016-10clusters');
P=cpar('NYC2018');
P.Operations.maxwait=Inf;
P.m=2000;
P.enlayeralg='no';
P.tripday=1;
P.modechoice=true;
P.e=1;

Pricing.tp=10;
Pricing.relocationcost=0.1;
Pricing.basetariff=0.25;
Pricing.alternative=0.25;
Pricing.VOT=15; % value of time
Pricing.pricingwaiting=0;
P.Pricing=Pricing;

P.Pricing.dynamic=0;
Res1=generalC(P,1,2)

P.Pricing.dynamic=1;
Res2=generalC(P,2,2)


%%

[   Res1.Sim.revenues;
    Res2.Sim.revenues]

[   Res1.Sim.relocationcosts;
    Res2.Sim.relocationcosts]

[   Res1.Sim.revenues-Res1.Sim.relocationcosts; 
    Res2.Sim.revenues-Res2.Sim.relocationcosts]

(Res2.Sim.revenues-Res2.Sim.relocationcosts)/((Res1.Sim.revenues-Res1.Sim.relocationcosts))

[   sum(Res1.Sim.chosenmode)/length(Res1.Sim.chosenmode);
    sum(Res2.Sim.chosenmode)/length(Res2.Sim.chosenmode)]

figure
plot(0:0.01:0.5,histc(Res2.Sim.offeredprices,0:0.01:0.5))

figure
x=0:0.2:10;
plot(x,histc(Res2.Sim.offeredprices,x))

%% sensitivity

Pmat=changeStruct(P,'Pricing.horizon',15:5:35,'Pricing.maxiter',1:6);
parfor i=1:length(Pmat)
    generalC(Pmat{i},1,-i);
end

% [~,~]=generateplotline3(mapscenario,[],varargin)

for i=1:length(Pmat)
    
    Res2=generalC(Pmat{i},1,2);
%     
% [   Res1.Sim.revenues;
%     Res2.Sim.revenues]
% 
% [   Res1.Sim.relocationcosts;
%     Res2.Sim.relocationcosts]
% 
% [   Res1.Sim.revenues-Res1.Sim.relocationcosts; 
%     Res2.Sim.revenues-Res2.Sim.relocationcosts]

    risulta(i)=(Res2.Sim.revenues-Res2.Sim.relocationcosts)/((Res1.Sim.revenues-Res1.Sim.relocationcosts));

% [   sum(Res1.Sim.chosenmode)/length(Res1.Sim.chosenmode);
%     sum(Res2.Sim.chosenmode)/length(Res2.Sim.chosenmode)]
    
    
end

risu=((reshape(risulta,6,5)-1)*100);

%%

% P=cpar('NYC2018');
% P.Operations.maxwait=Inf;
% P.m=2500;
% P.enlayeralg='no';
% P.tripday=1;
% P.modechoice=true;
% P.e=1;

Pricing.tp=10;
Pricing.relocationcost=0.1;
Pricing.basetariff=0.25;
Pricing.alternative=0.25;
Pricing.VOT=15; % value of time
Pricing.pricingwaiting=0;
Pricing.maxiter=5;
Pricing.horizon=20;

Pricing1=Pricing;
Pricing1.dynamic=0;
Pricing2=Pricing;
Pricing2.dynamic=1;
% Pricing1w=Pricing1;
% Pricing1w.pricingwaiting=1;
% Pricing2w=Pricing2;
% Pricing2w.pricingwaiting=2;

[P1,R1]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2500,'modechoice',true,'enlayeralg',"no",'e',1,'Pricing',Pricing1);
[P2,R2]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2500,'modechoice',true,'enlayeralg',"no",'e',1,'Pricing',Pricing2);

[P1,R1]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2000,'modechoice',true,'enlayeralg',"no",'e',1,'Pricing',Pricing1);
[P2,R2]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2000,'modechoice',true,'enlayeralg',"no",'e',1,'Pricing',Pricing2);

[P1,R1]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2000,'modechoice',true,'enlayeralg',"no",'e',1,'Pricing',Pricing1,'Pricing.pricingwaiting',1);
[P2,R2]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2000,'modechoice',true,'enlayeralg',"no",'e',1,'Pricing',Pricing2,'Pricing.pricingwaiting',1);

% [P1,R1]=generateplotline3('NYC2018',1,'tripfolder',"NYC2018_10wed_small_nodes",'tripday',1:10,'Operations.maxwait',Inf,'m',200,'modechoice',true,'enlayeralg',"no",'e',1,'Pricing',Pricing1);
% [P2,R2]=generateplotline3('NYC2018',1,'tripfolder',"NYC2018_10wed_small_nodes",'tripday',1:10,'Operations.maxwait',Inf,'m',200,'modechoice',true,'enlayeralg',"no",'e',1,'Pricing',Pricing2);


S1=[R1.Sim];
S2=[R2.Sim];

U1=sum([S1.revenues])-sum([S1.relocationcosts])
U2=sum([S2.revenues])-sum([S2.relocationcosts])

U2/U1

[mean([S1.revenues]) mean([S2.revenues])
mean([S1.relocationcosts]) mean([S2.relocationcosts])
mean([S1.revenues]-[S1.relocationcosts]) mean([S2.revenues]-[S2.relocationcosts])]

[mean(vertcat(S1.waiting)) mean(vertcat(S2.waiting))]

mean([S2.relocationcosts])/mean([S1.relocationcosts])

CPUt=[mean([R1.cputime]) mean([R2.cputime])]
(mean([R2.cputime])-mean([R1.cputime]))/24

% % share of people waiting more than 1 hour
% sum(vertcat(S1.waiting)>60)/length(vertcat(S1.waiting))
% sum(vertcat(S2.waiting)>60)/length(vertcat(S1.waiting))

boxplot([[S1.revenues];[S1.relocationcosts];[S2.revenues];[S2.relocationcosts]]')
boxplot([[S1.revenues];[S1.revenues]-[S1.relocationcosts];[S2.revenues];[S2.revenues]-[S2.relocationcosts]]')



%% waiting

DataFolder=setDataFolder();
k=1;
[A,Atimes,AbuckC,Distances]=loadTrips(P1{k});


wait1=accumarray(Atimes(:,1),full(R1(k).Sim.waiting))./histc(Atimes(logical(R1(k).Sim.chosenmode),1),1:max(Atimes(logical(R1(k).Sim.chosenmode),1)));
wait2=accumarray(Atimes(:,1),full(R2(k).Sim.waiting))./histc(Atimes(logical(R2(k).Sim.chosenmode),1),1:max(Atimes(logical(R2(k).Sim.chosenmode),1)));

x=linspace(0,24,length(wait1));
figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot(x,wait1)
plot(x,wait2)
xlim([0,24])
xlabel('hour')
ylabel('average wait time (minutes)')
legend({'constant';'dynamic'},'Location','NorthWest')
box on
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/Pricing/waittimes'],'-depsc2');

%%

DataFolder=setDataFolder();
figure('Units','centimeters','Position',[10,7,10,7])
hold on

h1=subplot(1,3,1)
boxplot([[S1.revenues];[S2.revenues]]')
% ylabel('$');
xticklabels({'C','D'})
title('               a')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');

h2=subplot(1,3,2)
boxplot([[S1.relocationcosts];[S2.relocationcosts]]')
xticklabels({'C','D'})
title('               b')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');

h3=subplot(1,3,3)
boxplot([[S1.revenues]-[S1.relocationcosts];[S2.revenues]-[S2.relocationcosts]]')
xticklabels({'C','D'})
title('               c')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');

p1=get(h1,'position');
p2=get(h2,'position');
p3=get(h3,'position');
height=p1(2)+p1(4)-p3(2);
h3=axes('position',[p1(1) p1(2) p2(3) height],'visible','off');
h_label=ylabel('$','visible','on');

set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/Pricing/results'],'-depsc2');













%%

Res2=generalC(Pmat{12},1,2);

% Res=R2(k);
% P0=P2{k};
Res=Res2;
P0=P;

DataFolder=setDataFolder();
k=1;
load([DataFolder 'scenarios/' P0.scenario],'T')
g=@(a,s) exp(s)./(exp(s)+a);        % value at s
d=@(a,s) a*exp(s)./((a+exp(s)).^2); % derivative at s
R=@(p,c) (exp(-p.*c)./(exp(-p.*c)+exp(-0.25*c))).*(p-0.1).*c; % net revenues at certain price and distance


figure('Units','centimeters','Position',[10,7,10,7])
hold on

h=4;
Prices=Res.Sim.prices(:,:,h*60/10);
changedprices=logical((Prices(:)~=0.25).*(Prices(:)~=0));
scatter(T(changedprices),Prices(changedprices),'.')

h=9;
Prices=Res.Sim.prices(:,:,h*60/10);
changedprices=logical((Prices(:)~=0.25).*(Prices(:)~=0));
scatter(T(changedprices)+60,Prices(changedprices),'.')

h=19;
Prices=Res.Sim.prices(:,:,h*60/10);
changedprices=logical((Prices(:)~=0.25).*(Prices(:)~=0));
scatter(T(changedprices)+120,Prices(changedprices),'.')

% theoretical optimal price points
for i=1:60
    optiprices(i)=bestp(i,0.1,0.25);
end
x=4:60;
p1=plot(x,optiprices(x),'k-');
plot(x+60,optiprices(x),'k-');
plot(x+120,optiprices(x),'k-');

xlim([0,180])
xticks([0:20:180])
xticklabels(num2str([0,20,40,0,20,40,0,20,40,60]'))
ylabel('price per minute ($)')
xlabel('OD pair distance (min)')
text(20,0.5,'4:00 am')
text(20+50,0.5,'9:00 am')
text(20+100,0.5,'7:00 pm')
box on
grid on
legend(p1,{'theoretical best price'})
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
% print([DataFolder 'figures/Pricing/pricevariation'],'-depsc2');
