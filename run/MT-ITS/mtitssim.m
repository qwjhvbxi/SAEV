


P=cpar('NYC2018');
P.Operations.maxwait=Inf;
P.m=5000;
P.enlayeralg='no';
P.tripday=3;
P.modechoice=true;
P.e=1;

Pricing.tp=30;
Pricing.relocationcost=0.1;
Pricing.basetariff=0.25;
Pricing.alternative=0.25;
Pricing.VOT=15; % value of time
Pricing.pricingwaiting=0;

P.Pricing=Pricing;

P.Pricing.dynamic=0;
Res1=generalC(P,1,2)

P.Pricing.dynamic=1;
Res2=generalC(P,1,2)

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

%%

addpath plots functions utilities
% [P1,R1]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',5000,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',false,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P2,R2]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',5000,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25,'TransportLayer.alternative',0.25);
% [P1,R1]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2500,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',false,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P2,R2]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2500,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P1,R1]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2500,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',false,'pricingwaiting',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P2,R2]=generateplotline3('NYC2018',1,'tripday',1:10,'Operations.maxwait',Inf,'m',2500,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',true,'pricingwaiting',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P1,R1]=generateplotline3('NYC2018',1,'tripfolder',"NYC2018_10wed_small",'tripday',1:10,'Operations.maxwait',Inf,'m',500,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',false,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P2,R2]=generateplotline3('NYC2018',1,'tripfolder',"NYC2018_10wed_small",'tripday',1:10,'Operations.maxwait',Inf,'m',500,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P1,R1]=generateplotline3('NYC2018',1,'tripfolder',"NYC2018_10wed_small",'tripday',1:10,'Operations.maxwait',Inf,'m',250,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',false,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P2,R2]=generateplotline3('NYC2018',1,'tripfolder',"NYC2018_10wed_small",'tripday',1:10,'Operations.maxwait',Inf,'m',250,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P1,R1]=generateplotline3('NYC2018',1,'tripfolder',"NYC2018_10wed_small",'tripday',1:10,'Operations.maxwait',Inf,'m',250,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',false,'pricingwaiting',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
% [P2,R2]=generateplotline3('NYC2018',1,'tripfolder',"NYC2018_10wed_small",'tripday',1:10,'Operations.maxwait',Inf,'m',250,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',true,'pricingwaiting',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);

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

%%

Pricing.relocationcost=0.1;
Pricing.basetariff=0.25;
Pricing.alternative=0.25;
Pricing.tp=30;
Pricing.pricingwaiting=false;
Pricing.dynamic=true;
Pricing.VOT=15;

addpath plots functions utilities
[P1,R1]=generateplotline3('NYC2018',-1,'tripday',1,'Operations.maxwait',Inf,'m',5000,'modechoice',true,'Pricing',Pricing,'enlayeralg',"no",'Pricing.dynamic',false);
[P2,R2]=generateplotline3('NYC2018',-1,'tripday',1,'Operations.maxwait',Inf,'m',5000,'modechoice',true,'Pricing',Pricing,'enlayeralg',"no",'Pricing.dynamic',true);


%%

DataFolder=setDataFolder();
k=1;
[A,Atimes,AbuckC,Distances]=generateGPStrips(P1{k});


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

set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/TRB/waittimes'],'-depsc2');

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
print([DataFolder 'figures/TRB/results'],'-depsc2');


%% probabilities

DataFolder=setDataFolder();
figure('Units','centimeters','Position',[10,7,10,7])
hold on
p=0:0.01:0.5;
D=[5 10 20];
cc=lines(3);
for k=1:3
    d=D(k);
    yyaxis right
    J1(k)=plot(p,f2(p,d),'--','Color',cc(k,:),'LineWidth',1);
    yyaxis left
    J2(k)=plot(p,f2(p,d).*p,':','Color',cc(k,:),'LineWidth',1);
    J3(k)=plot(p,f2(p,d).*p-f2(p,d)*0.1,'-','LineWidth',1.5,'Color',cc(k,:));
end
xlabel('price per minute')
ylabel('$ per minute')
yyaxis right
ylabel('probability of acceptance')
text(0.4,0.82,'d=5')
text(0.43,0.62,'d=10')
text(0.4,0.4,'d=20')
legend([J1(1),J2(1),J3(1)],{'probability','revenue','profits'})
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/TRB/revenue'],'-depsc2');


%% linearization intervals

m.gamma_p=0.25;
m.gamma_r=0.1;

g=@(a,s) exp(s)./(exp(s)+a);        % value at s
f2=@(s,d) exp(-s.*d)./(exp(-s.*d)+exp(-m.gamma_p*d)); % demand acceptance probability as a function of price and distance
q=@(d,x) -(log(x.*exp(-m.gamma_p*d)./(1-x)))./d;
Points=g(exp(0),-3.5:3.5);          % probability linearization intervals (7 intervals, 8 limits)
D=[5 10 20];
p=0:0.01:0.5;
cc=lines(3);

DataFolder=setDataFolder();
figure('Units','centimeters','Position',[10,7,10,7])
hold on
for k=1:length(D)
    plot(p,f2(p,D(k)),'k--')
    w=0; % limits for p
    for w=-3:3        
        % probabilities of trip acceptance at the limits
        fmin=1-Points(5+w);
        fmax=1-Points(4+w);
        % limits to price
        pmin=q(D(k),fmax);
        pmax=q(D(k),fmin);
        
        line([pmin,pmax],[fmax,fmin],'LineWidth',1,'Color',cc(k,:));
        scatter([pmin,pmax],[fmax,fmin],'ko');%,'MarkerColor',cc(k,:));
    end
end
xlim([0,0.5])
xlabel('price per minute')
ylabel('probability of acceptance')
text(0.025,0.7,'\it d=5','FontName','Times')
text(0.05,0.85,'\it d=10','FontName','Times')
text(0.2,0.82,'\it d=20','FontName','Times')

w=1;
fmin=1-Points(5+w);
fmax=1-Points(4+w);
pmin=q(D(k),fmax);
pmax=q(D(k),fmin);
line([[0;pmin],[0;pmax]],[[fmax;fmax],[fmin;fmin]],'LineStyle',':','Color','k')
line([[pmin;pmin],[pmax;pmax]],[[0;fmax],[0;fmin]],'LineStyle',':','Color','k')
text(0.01,fmin,'{\it p^{min}}','FontName','Times')
text(0.01,fmax,'{\it p^{max}}','FontName','Times')
text(pmin-0.03,0.05,'{\it z^{min}}','FontName','Times')
text(pmax-0.01,0.05,'{\it z^{max}}','FontName','Times')

set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/TRB/linearization2'],'-depsc2');


%% prices graph

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

h=8;
Prices=Res.Sim.prices(:,:,h);
changedprices=(Prices(:)~=0.25);
scatter(T(changedprices),Prices(changedprices),'.')

h=3;
Prices=Res.Sim.prices(:,:,h);
changedprices=(Prices(:)~=0.25);
scatter(T(changedprices),Prices(changedprices),'.')

% theoretical optimal price points
for i=1:50
    optiprices(i)=bestp(i);
end
x=4:50;
plot(x,optiprices(x),'k-');

xlim([0,50])
ylabel('price per minute ($)')
xlabel('OD pair distance (min)')
text(0.025,0.7,'d=5')
text(0.05,0.85,'d=10')
text(0.2,0.82,'d=20')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
legend({'7:00-8:00 am';'2:00-3:00 am';'theoretical best price'})
% print([DataFolder 'figures/TRB/pricevariation2'],'-depsc2');



