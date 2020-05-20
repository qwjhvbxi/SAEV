%% simulations for JSER paper



%% simulations

% addpath plots
% [P1,R1]=generateplotline3('Tokyo189',[],'Operations.maxwait',Inf,'gridday',1:5);

P=cpar('Tokyo189');
% P.Operations.maxwait=Inf;
for k=1:365
    
    P.gridday=k;
    R(k)=generalC(P,1,-k);
    P.Operations.initialsoc=R(k).Sim.q(end,:);
    uinit=max(1,min(189,R(k).Sim.u(end,:)+sum(R(k).Internals.v(721:end,:))+sum(R(k).Internals.w(721:end,:))));
    P.Operations.uinit=uinit;
    
end

q=[];
for k=1:365
    
    q=[q;mean(R(k).Sim.q,2)];
%     
%     R(k).Sim.u=uint8(R(k).Sim.u); % final destination of vehicles (station) [tsim x m]
%     R(k).Sim.q=single(R(k).Sim.q); % state of charge 
%     R(k).Sim.e=sparse(R(k).Sim.e);
%     R(k).Internals.v=sparse(R(k).Internals.v);
%     R(k).Internals.w=sparse(R(k).Internals.w);
    
end
plot(q)


% calculate baseline cost
load('data/eleprices/TokyoDA-FY2017-Reduced.mat','x');
P0=cpar('Tokyo189');
P0.enlayeralg='no';
P0.Operations.initialsoc=1;
R0=generalC(P0,1,2);

cost0=(sum(R0.Sim.e/60*P0.e,2)')*repelem(x,P0.beta,1);


costComp=[cost0',[R(:).cost]'];

DataFolder=setDataFolder();

figure
plot(costComp)

% Cost comparison
figure('Units','centimeters','Position',[10,7,10,7])
histVec=-6*10^6:10^5:1.2*10^6;
hold on
plot(histVec(1:end-1),histcounts(costComp(:,1),histVec))
plot(histVec(1:end-1),histcounts(costComp(:,2),histVec))
xlabel('daily cost (yen)')
ylabel('number of days')
legend({'unscheduled','scheduled'},'Orientation','vertical','Location','NorthWest')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/JSER/CostComparison'],'-dpng','-r300');


% charging power and electricity price

z=linspace(0,24,721);

figure('Units','centimeters','Position',[10,7,10,7])
Day=311;% 7 Nov 2017
hold on
yyaxis right
plot(z(1:end-1),R(Day).Params.elep(1:720))
ylabel('electricity price (yen/kWh)')
yyaxis left
stairs(linspace(0,23.5,48),((R(Day).Internals.zmacro(1,1:48)-R(Day).Internals.zmacro(2,1:48)).*R(Day).Internals.zmacro(3,1:48))'*60/30/1000)
plot(z(1:end-1),sum(R(Day).Sim.e,2)/1000,'k-')
xlim([0,24])
xlabel('hour')
ylabel('power (MW)')
lgnd=legend({'EL power','TL power','price'},'Orientation','vertical','Location','NorthWest');
lgnd.BoxFace.ColorType='truecoloralpha';
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
% print([DataFolder 'figures/JSER/charging'],'-dpng','-r300');
print('charging','-dpng','-r300');


figure('Units','centimeters','Position',[10,7,10,7])
plot(z(1:end-1),sum(R0.Sim.e,2))
xlabel('daily cost (yen)')
ylabel('number of days')
legend({'unscheduled','scheduled'},'Orientation','vertical','Location','NorthWest')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/JSER/charging'],'-dpng','-r300');

sum(costComp)



return


%% sensitivity analysis

addpath plots
[P1,R1]=generateplotline3('Tokyo189',[],'Operations.maxwait',Inf,'m',[1500 2000 2500 3000 3500 4000 4500 5000]);


% joy division graph for relocation
figure
hold on
t=5;
for j=1:length(R1)
    plot(sum(reshape(R1(length(R1)+1-j).Sim.relodist+j*200,t,720/t)))
end


% joy division graph for SOC
figure
hold on
t=5;
for j=1:length(R1)
    plot(mean(reshape(mean(R1(length(R1)+1-j).Sim.q(1:720,:),2)+j*0.1,t,720/t)))
end




return

%% generate Tokyo scenario 

DataFolder=setDataFolder();

load([DataFolder 'input_files/Tokyo/sources/tripsCell.mat'],'Trips');
load([DataFolder 'input_files/Tokyo/sources/aev18data.mat'],'areas','chargingstnodes','gridcoord','tph');

% each trip in Trips is generated according to a gaussian centered in the
% hour. Weights are between 16 and 210. I divide by 15 and round, so I get
% about 188,000 trips

%% create A, Atimes

ScaleFactor=15;
A=[];
Atimes=[];
for h=1:24
    
    TripWeights=max(1,round(Trips{h}(:,1)./ScaleFactor));
    HourlyTrips=repelem(Trips{h}(:,2:3),TripWeights,1);
    RandMinutes=round((h-1)*60+30+randn(length(HourlyTrips),1)*30);
    A=[A;HourlyTrips];
    Atimes=[Atimes;[RandMinutes RandMinutes]];
    
end

% bring trips out of either side of the day on the other side
Atimes(Atimes>1440)=Atimes(Atimes>1440)-1440;
Atimes(Atimes<1)=Atimes(Atimes<1)+1440;

% remove trips starting and ending in same node
SameNode=A(:,1)==A(:,2);
A(SameNode,:)=[];
Atimes(SameNode,:)=[];

% save([DataFolder 'trips/Tokyo2008_1day_189k.mat'],'A','Atimes'); % with same node trips
save([DataFolder 'trips/Tokyo2008_1day_162k.mat'],'A','Atimes');


%% create scenario with less stations

limits=[10 30 15 35];
selection=logical((gridcoord(:,1)>limits(1)).*(gridcoord(:,1)<limits(2)).*(gridcoord(:,2)>limits(3)).*(gridcoord(:,2)<limits(4)));
AllowedStations=find(selection);

% create vector to assign new relative indeces to nodes
NewIndeces=zeros(514,1);
NewIndeces(AllowedStations)=1:length(AllowedStations);

% figure
% scatter(gridcoord(:,1),gridcoord(:,2))
% figure
% scatter(gridcoord(selection,1),gridcoord(selection,2))

TripSelection=logical(ismember(A(:,1),AllowedStations).*ismember(A(:,2),AllowedStations));
A=NewIndeces(A(TripSelection,:));
Atimes=Atimes(TripSelection,:);
save([DataFolder 'trips/Tokyo2008_1day_48k.mat'],'A','Atimes');


%% create C, T (from Euclidean distance)

Eta=1.5;
T=max(1,round(Eta*sqrt((gridcoord(:,1)-gridcoord(:,1)').^2+(gridcoord(:,2)-gridcoord(:,2)').^2)));

C=gridcoord;
save('data/scenarios/Tokyo514.mat','C','T');

C=gridcoord(selection,:);
T=T(selection,selection);
save('data/scenarios/Tokyo189.mat','C','T')


%% map of Tokyo

load('data/scenarios/Tokyo189.mat');

TokyoMap=imread([DataFolder 'input_files/Tokyo/maps/TokyoStreetMap.png']);
limits=[10 30 15 35];
TokyoMapFlipped=flipdim(TokyoMap,1);

figure
hold on
image(0:40,0:40,TokyoMapFlipped);
scatter(C(:,1),C(:,2),'ko','MarkerFaceColor','k')
% xticks(limits(1):5:limits(2))
% xticklabels(0:5:20)
% yticks(limits(3):5:limits(4))
% yticklabels(0:5:20)
axis equal off
axis(limits)
print([DataFolder 'figures/JSER/TokyoMap189'],'-dpng','-r300');





