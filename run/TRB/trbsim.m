

% P=cpar('Tokyo189');
P=cpar('NYC2018');
% P=cpar('NYC2016');
P.Operations.maxwait=Inf;
P.m=2500;
P.TransportLayer.tp=30;
P.enlayeralg='no';
P.TransportLayer.relocationcost=0.1;
P.TransportLayer.basetariff=0.25;
P.tripday=3;

P.pricing=false;
Res1=generalC(P,1,2)

P.pricing=true;
Res2=generalC(P,1,2)

[   Res1.Sim.revenues;
    Res2.Sim.revenues]

[   Res1.Sim.relocationcosts;
    Res2.Sim.relocationcosts]

[Res1.Sim.revenues-Res1.Sim.relocationcosts; Res2.Sim.revenues-Res2.Sim.relocationcosts]
(Res2.Sim.revenues-Res2.Sim.relocationcosts)/((Res1.Sim.revenues-Res1.Sim.relocationcosts))

[sum(Res1.Sim.chosenmode)/length(Res1.Sim.chosenmode);
sum(Res2.Sim.chosenmode)/length(Res2.Sim.chosenmode)]

figure
plot(0:0.01:0.5,histc(Res2.Sim.offeredprices,0:0.01:0.5))

%%

addpath plots
[P1,R1]=generateplotline3('NYC2018',2,'tripday',1:10,'Operations.maxwait',Inf,'m',5000,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',false,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
[P2,R2]=generateplotline3('NYC2018',2,'tripday',1:10,'Operations.maxwait',Inf,'m',5000,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
[P1,R1]=generateplotline3('NYC2018',2,'tripday',1:10,'Operations.maxwait',Inf,'m',2500,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',false,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);
[P2,R2]=generateplotline3('NYC2018',2,'tripday',1:10,'Operations.maxwait',Inf,'m',2500,'TransportLayer.tp',30,'enlayeralg',"no",'pricing',true,'TransportLayer.relocationcost',0.1,'TransportLayer.basetariff',0.25);

S1=[R1.Sim];
S2=[R2.Sim];

U1=sum([S1.revenues])-sum([S1.relocationcosts])
U2=sum([S2.revenues])-sum([S2.relocationcosts])

U2/U1

[mean([S1.revenues]) mean([S2.revenues])
mean([S1.relocationcosts]) mean([S2.relocationcosts])
mean([S1.revenues]-[S1.relocationcosts]) mean([S2.revenues]-[S2.relocationcosts])]

mean([R1.cputime])
mean([R2.cputime])
(mean([R2.cputime])-mean([R1.cputime]))/24

% share of people waiting more than 1 minute
sum(vertcat(S1.waiting)>60)/length(vertcat(S1.waiting))
sum(vertcat(S2.waiting)>60)/length(vertcat(S1.waiting))

boxplot([[S1.revenues];[S1.relocationcosts];[S2.revenues];[S2.relocationcosts]]')
boxplot([[S1.revenues];[S1.revenues]-[S1.relocationcosts];[S2.revenues];[S2.revenues]-[S2.relocationcosts]]')

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
text(0.4,0.82,'c=5')
text(0.43,0.62,'c=10')
text(0.4,0.4,'c=20')
legend([J1(1),J2(1),J3(1)],{'probability','revenue','profits'})
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/TRB/revenue'],'-depsc2');


%% linearization intervals

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
    u=0.5;
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
text(0.025,0.7,'c=5')
text(0.05,0.85,'c=10')
text(0.2,0.82,'c=20')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/TRB/linearization'],'-depsc2');







