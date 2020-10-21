%% comparison between optimal and simplified

% gridday=13;
% 
% P1=cpar('NYC2016-small','opti');
% P1.gridfile='NY_DA_2016';
% P1.gridday=gridday;
% P1.Operations.initialsoc=0.5;
% Res1=generalC(P1,1,2)
% 
% P2=cpar('NYC2016-small','simplified');
% P2.gridfile='NY_DA_2016';
% P2.gridday=gridday;
% P2.Operations.maxwait=Inf;
% P2.Operations.initialsoc=0.5;
% % P2.Tech.cyclingcost=100/4000;
% Res2=generalC(P2,1,2)



% %% 1 week
% 
% Period=25:31;
% gridoffset=0; % same year
% 
% P=cpar('NYC2016-small','opti');
% P.tripfolder='NYC2016-small'; % trips in 2019
% P.Operations.maxwait=Inf;
% P.gridfile='NY_DA_2016';
% [S1,R1]=multiDaySim(Period,P,gridoffset);
% 
% P=cpar('NYC2016-small','simplified');
% P.tripfolder='NYC2016-small'; % trips in 2019
% P.Operations.maxwait=Inf;
% P.gridfile='NY_DA_2016';
% [S2,R2]=multiDaySim(Period,P,gridoffset);


%% 1 week

addpath functions utilities plots
DataFolder=setDataFolder();
% Days=18:24;
Days=25:31;
m=30;
scenario='NYC2016-small';

% optimal
Transport.thor=8;          % horizon (time steps)
Transport.rho1=0.01;        % weight of secondary objective
Transport.rho2=0.01;        % weight of charging objective for electricity price
Transport.rho3=0.01;        % weight of charging objective for SOC
Transport.rho4=0.000001;    % weight for fixed charge

Trips={'NYC2016-small_200','NYC2016-small_150','NYC2016-small_100','NYC2016-small_75','NYC2016-small_50'};


% load or create start variables to avoid saving files in parallel
% computing
load([DataFolder 'scenarios/' scenario],'T');
n=size(T,1);
StartFile=[DataFolder 'temp/StartFile_' num2str(n) '-' num2str(m) '.mat'];
if exist(StartFile,'file')
    load(StartFile,'StartSoc','StartPos');
else
    StartSoc=ones(1,P.m)*0.7;
    StartPos=randi(n,1,P.m);
    save(StartFile,'StartSoc','StartPos');
end

% Transport2.tx=10;  % minutes
% Transport2.ts=24; 
% Transport2.tr=20; 
% Transport2.bmin=0;

% create P cell
K=length(Trips);
for k=1:K
    
    % common parameters
    Pmat{k}=cpar(scenario);
    Pmat{k}.m=m;
    Pmat{k}.Operations.maxwait=Inf;
    Pmat{k}.trlayeralg='opti';
    Pmat{k}.TransportLayer=Transport;
    Pmat{k}.gridfile='NY_DA_2016';
    Pmat{k}.tripfolder=Trips{k};
    
    
    % simplified
    Pmat{K+k}=cpar(scenario);
    Pmat{K+k}.m=m;
    Pmat{K+k}.e=2;
    Pmat{K+k}.Operations.maxwait=Inf;
    Pmat{K+k}.trlayeralg='simplified';
    % Pmat{K+k}.TransportLayer=Transport2;
    Pmat{K+k}.gridfile='NY_DA_2016';
    Pmat{K+k}.tripfolder=Trips{k};
    
end

gridoffset=0;

parfor k=1:length(Pmat)
    multiDaySim(Days,Pmat{k},gridoffset)
end


for k=1:length(Pmat)
    [S1,R1]=multiDaySim(Days,Pmat{k},gridoffset)
    S(k)=S1;
%     R(k)=R1;
end


%% tests

[S1,R1]=multiDaySim(Days(1),Pmat{4},gridoffset,1)
[S2,R2]=multiDaySim(Days(1),Pmat{8},gridoffset,1)

figure
hold on
plot(sum(R1.Sim.e,2))
plot(sum(R2.Sim.e,2))
plotta(R1,'power')
plotta(R2,'power')

%% plot cputime

linevect={'x-','s-','v-','^-','d-','>--','<--','o--','x--'};
figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot(x,sum([S(1:4).cputime]),linevect{2})
plot(x,sum([S(5:8).cputime]),'^--')
% semilogy(x,[sum([S(1:4).cputime])'  sum([S(5:8).cputime])'],'x-');
set(gca, 'YScale', 'log')
grid on
ylabel('CPU time (s)')
xticks(x)
xlabel('total trips (1 week)')
lgnd=legend({'exact','aggregated'},'Orientation','vertical','Location','East');
% lgnd.BoxFace.ColorType='truecoloralpha';
% lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/Energy/comparison_cputime'],'-depsc2');

%% plot wait cost

linevect={'x-','s-','v-','^-','d-','>--','<--','o--','x--'};
N=cellfun('length',{S.waiting});
x=N(1:4);

z=linspace(0,24,721);

figure('Units','centimeters','Position',[10,7,10,7])
hold on
L1=plot(x,max([S(1:4).cost]),linevect{2});
L2=plot(x,max([S(5:8).cost]),linevect{7});
ylabel('charging cost ($)')

yyaxis right
L3=plot(x,max([S(1:4).peakwait]),linevect{1});
L4=plot(x,max([S(5:8).peakwait]),linevect{6});
% plot(x,max([S(1:4).avgwait]),linevect{1});
% plot(x,max([S(5:8).avgwait]),linevect{6});
ylabel('peak wait time (min)')
ylim([0,13])


xticks(x)
xlabel('total trips (1 week)')
lgnd=legend([L1,L2,L3,L4],{'cost (exact)','cost (aggregated)','peak wait (exact)','peak wait (aggregated)'},'Orientation','vertical','Location','SouthEast');
lgnd.BoxFace.ColorType='truecoloralpha';
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/Energy/comparison_costpeak'],'-depsc2');

%%

% figure
% plot(reshape([S.totavgwait],4,2))
% figure
% semilogy(reshape([S.totavgwait],4,2))
% 
% figure
% plot([S(1:4).soc])
% 
% figure
% plot([S(5:8).soc])

% [Pz,Rz]=generateplotline3(scenario,[],'Operations.maxwait',Inf,'gridfile',"NY_DA_2016",'tripfolder',"NYC2016-small_50",'m',30:5:120,'tripday,gridday',30)

%%

if 0

    %% 5 wednesdays 

    Days=13:7:13+7*4;
    % Days=13

    % optimal
    Transport.thor=8;          % horizon (time steps)
    Transport.rho1=0.01;        % weight of secondary objective
    Transport.rho2=0.01;        % weight of charging objective for electricity price
    Transport.rho3=0.01;        % weight of charging objective for SOC
    Transport.rho4=0.000001;    % weight for fixed charge

    % optimal
    [P1,R1]=generateplotline3('NYC2016-small',[],'Operations.maxwait',Inf,'trlayeralg',"opti",'TransportLayer',Transport,'gridfile',"NY_DA_2016",'tripfolder',["NYC2016-small_150","NYC2016-small_100"],'tripday,gridday',Days);

    % simplified
    [P2,R2]=generateplotline3('NYC2016-small',[],'Operations.maxwait',Inf,'trlayeralg',"simplified",'gridfile',"NY_DA_2016",'tripfolder',"NYC2016-small_150",'tripday,gridday',Days);

    % cost comparison
    sum([R1.cost])
    sum([R2.cost])

    % waiting comparison
    Sims2=[R2.Sim];
    Totwaiting=full(vertcat(Sims2.waiting));

    sum(Totwaiting>0)/length(Totwaiting)
    mean(Totwaiting(Totwaiting>0))

    sum([R1.cost])
    sum([R2.cost])

end


if 0

    d=1;
    
    Res1=R1(d);
    Res2=R2(d);
    Par2=P2{d};
    
    %% plots

    z=linspace(0,24,721);

    figure('Units','centimeters','Position',[10,7,10,7])
    hold on

    yyaxis right
    plot(z(1:end-1),Res1.Params.elep(1:720))
    ylabel('electricity price ($/MWh)')

    yyaxis left
    plot(z(1:end-1),sum(Res1.Sim.e,2)/1000)
    plot(z(1:end-1),sum(Res2.Sim.e,2)/1000,'k-')
    xlim([0,24])
    xlabel('hour')
    ylabel('power (MW)')
    lgnd=legend({'exact','aggregated','price'},'Orientation','vertical','Location','North');
%     lgnd.BoxFace.ColorType='truecoloralpha';
%     lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
    set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
    print([DataFolder 'figures/Energy/comparison'],'-depsc2');
    
    
    
    figure('Units','centimeters','Position',[10,7,10,7])
    hold on
    plot(z(1:end),mean(Res1.Sim.q,2))
    plot(z(1:end),mean(Res2.Sim.q,2),'k-')
    xlim([0,24])
    xlabel('hour')
    ylabel('SOC')
    lgnd=legend({'exact','aggregated'},'Orientation','vertical','Location','South');
%     lgnd.BoxFace.ColorType='truecoloralpha';
%     lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
    set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
    print([DataFolder 'figures/Energy/comparison_soc'],'-depsc2');

    
    Reqs=zeros(1440,length(R1));
    Waits=zeros(1440,length(R1));

    for d=1:length(R1)

        Res1=R1(d);
        Res2=R2(d);
        Par2=P2{d};

    %     load([DataFolder 'trips/' Par2.tripfile])
        load([DataFolder 'trips/' Par2.tripfolder '/d' num2str(Par2.tripday) '.mat'])

        Atimes(Atimes(:,1)==0,1)=1;
        Reqs(:,d)=histc(Atimes(:,1),1:1440);
        Waits(:,d)=accumarray([Atimes(:,1);1440],[Res2.Sim.waiting;0]);

    end
    
%     plot(Waits(:,d))
%     plot(Reqs(:,d))
% 
%     Waits10min=sum(reshape(Waits(:,d),10,144));
%     Reqs10min=sum(reshape(Reqs(:,d),10,144));


    Bucket=60;
    Waits10min=sum(reshape(sum(Waits,2),Bucket,1440/Bucket));
    Reqs10min=sum(reshape(sum(Reqs,2),Bucket,1440/Bucket));

    x=linspace(0,24,1440/Bucket);

    figure('Units','centimeters','Position',[10,7,10,7])
    hold on
    line([0,24],[0,0])
    stairs(x,Waits10min./Reqs10min,'k:')
    xlim([0,24])
    xlabel('hour')
    ylabel('average wait time (min.)')
    lgnd=legend({'exact','aggregated'},'Orientation','vertical','Location','NorthWest');
    lgnd.BoxFace.ColorType='truecoloralpha';
    lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
    set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
    % print([DataFolder 'figures/JSER/charging'],'-dpng','-r300');


    %%

    Res=Res1;

    figure('Units','centimeters','Position',[10,7,10,7])
    hold on
    yyaxis right
    plot(z(1:end-1),Res.Params.elep(1:720))
    ylabel('electricity price (yen/kWh)')
    yyaxis left
    stairs(linspace(0,23.5,48),((Res.Internals.zmacro(1,1:48)-Res.Internals.zmacro(2,1:48)).*Res.Internals.zmacro(3,1:48))'*60/30/1000)
    plot(z(1:end-1),sum(Res.Sim.e,2)/1000,'k-')
    xlim([0,24])
    xlabel('hour')
    ylabel('power (MW)')
    lgnd=legend({'EL power','TL power','price'},'Orientation','vertical','Location','NorthWest');
    lgnd.BoxFace.ColorType='truecoloralpha';
    lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.6]');
    set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
    % print([DataFolder 'figures/JSER/charging'],'-dpng','-r300');

    figure('Units','centimeters','Position',[10,7,10,7])
    plot(z(1:end-1),sum(R0.Sim.e,2))
    xlabel('daily cost (yen)')
    ylabel('number of days')
    legend({'unscheduled','scheduled'},'Orientation','vertical','Location','NorthWest')
    set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
    % print([DataFolder 'figures/JSER/charging'],'-dpng','-r300');

end

% %% debug simulation
% 
% P1=cpar('NYC2016-small2','simplified');
% P1.Operations.maxwait=Inf;
% Res1=generalC(P1,2,2)
% P=cpar('NYC2016-small2','opti');
% Res2=generalC(P,2,2)
% 
% return