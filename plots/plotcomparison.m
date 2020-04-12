%% plot comparisons


P(1).tsim=1440/P(1).e;    % length of simulation (time steps) (one day)
P(1).mtsim=P(1).tsim/P(1).macrosteps; % length of macro simulation
P(1).speed=P(1).speedkmh/60*P(1).e;  % speed km/h -> (units/time step)
P(1).ac=round(P(1).chargekw/P(1).battery/60*P(1).e,3);      % charge rate
P(1).ad=P(1).consumption*P(1).speedkmh/P(1).battery/60*P(1).e;      % discharge rate (consumption*speed/batery/minutes)
P(1).T=P(1).tsim/P(1).macrosteps;
        

%% comparison between n groups (for example different P.types 1/6)

% manual selection
groups=[39:48
        29 , 49:52 , 34:38%29:38 % replaced simulations with mistaken same arrivalseed
        ];
Z=cell(size(groups,1),1)
for i=1:size(groups,1)
    Z{i}=outdesc(groups(i,:),0);
end

% automatic
groups=zeros(3,20);
Z=cell(size(groups,1),1);
Z{1}=outdesc([],[],'eleprofseed',eleplist,'type',1);
Z{2}=outdesc([],[],'eleprofseed',eleplist,'type',6);

%% variables
Z=vertcat(Z{:})
types=reshape([Z.type],size(Z))
totalcost=reshape([Z.totalcost],size(Z))
totalcostfinalsoc=reshape([Z.totalcostfinalsoc],size(Z))
totalwait=[Z.waittimes]
avgwait=reshape(totalwait(8,:),size(Z))
peakwait=reshape([Z.waitmovavg10minPeak],size(Z))%reshape(totalwait(9,:),size(Z))
peakwaitabsolute=reshape(totalwait(7,:),size(Z));
wait975=reshape(totalwait(6,:),size(Z)); % wait times for 97.5percentile
totalcomptime=[Z.comptime]
avgcomptime=reshape(totalcomptime(5,:),size(Z))


%% total wait graph

% moving average
e=2;
tlen=720;
waiting=vertcat(Z(1,:).waiting);
binsize=round(60/e);
halfbin=floor(binsize/2);
waitingprof4=zeros(tlen,1);
for k2=1:tlen
    
    % waitingprof4(k2)=mean(    waiting(   logical((waiting(:,1)>=k2) - (waiting(:,1)>k2+binsize))   ,3)   ,'omitnan');
    waitingprof4(k2)=mean(    waiting(   logical((waiting(:,1)>=k2-halfbin).*(waiting(:,1)<=k2+halfbin))   ,3)   ,'omitnan');
    
end
waitingprof4(isnan(waitingprof4))=0; % for bins without arrivals, assume 0 waiting
figure
plot(waitingprof4);


%% planned vs effective availability (second group)

groupid=1;

dktripeff=[Z(groupid,:).dktripeff];
dktrip=[Z(groupid,:).dktrip];
dktrip=dktrip(1:size(dktripeff,1),:);
totdk=P(1).macrosteps*P(1).m;

% square error
figure
plot(mean((dktrip-dktripeff).^2,2))

% MAPE
MAPE=mean((dktripeff-dktrip)./dktripeff,2,'omitnan')*100
deviationperc=(-dktripeff+dktrip)./totdk*100
mean(deviationperc,2,'omitnan')
figure
plot(MAPE)
figure
plot(deviationperc)

% available vehicles: confidence intervals 
x=deviationperc';
SEM = std(x)./sqrt(size(x,1));               % Standard Error
ts = tinv([0.05  0.95],size(x,1)-1);      % T-Score
CI = mean(x) + ts'*SEM;                      % Confidence Intervals

% avilable vehicles: errors and confidence intervals
figure('pos',[500 200 400 300])
plot(linspace(0,24,P(1).T),mean(x)' )
hold on
plot(linspace(0,24,P(1).T),CI','k--')
line([0,24],[0,0],'LineStyle',':','Color',[0.25,0.25,0.25])
xlim([0,24])
xticks((0:4:24))
xlabel('hour','FontSize',12)
ylabel('error %','FontSize',12)
% saveas(gcf,'figs/connectederror.eps','epsc')

% expected vs real distance traveled comparison
figure
plot(dktrip,'k-')
hold on
plot(dktripeff,'r--')

%% expected vs real distance traveled difference
figure
plot(linspace(0,24,P(1).T),max(0,P(1).m*P(1).macrosteps-[dktrip(1:P(1).T) dktripeff]))
xlim([0,24])
xticks((0:4:24))


%% waiting

figure
avgwaittimesteps=[Z.waitmovavg10min];
avgwaittimesteps(isnan(avgwaittimesteps))=0;
avgwaittimesteps=reshape(avgwaittimesteps,[length(avgwaittimesteps),size(Z')]);
avgwaittimesteps=squeeze(mean(avgwaittimesteps,2));
plot(avgwaittimesteps)



%% wait/costs plots

% average wait times
figure('pos',[500 200 400 300])
boxplot(avgwait')
ylim([0 max(max(avgwait))*1.1])
% xlabel('case','FontSize',12)
ylabel('minutes','FontSize',12)
xticklabels({'unscheduled';'scheduled'})
saveas(gcf,'figs/avgwait.eps','epsc')

% peak wait times
figure('pos',[500 200 400 300])
boxplot(peakwait')
ylim([0 max(max(peakwait))*1.1])
% xlabel('case','FontSize',12)
ylabel('minutes','FontSize',12)
xticklabels({'unscheduled';'scheduled'})
saveas(gcf,'figs/peakwait.eps','epsc')

% costs comparison
figure('pos',[500 200 400 300])
boxplot(totalcostfinalsoc'/1000)
hold on
weightedcosts=totalcostfinalsoc*elepprob;
plot(1:2,weightedcosts/1000,'xk')
offset1=0.2;
text(1+offset1,weightedcosts(1)/1000+offset1,sprintf('%0.0f yen',weightedcosts(1)))
text(2-3*offset1,weightedcosts(2)/1000+offset1,sprintf('%0.0f yen',weightedcosts(2)))
% annotation('textbox',1,weightedcosts(1),sprintf('%0.0f',weightedcosts(1)))
% annotation('textbox',2,weightedcosts(2),'String',sprintf('%0.0f',weightedcosts(2)))
ylim([0 max(max(totalcostfinalsoc))/1000*1.1])
% xlabel('case','FontSize',12)
ylabel('thousand yen','FontSize',12)
xticklabels({'unscheduled';'scheduled'})
saveas(gcf,'figs/totalcostfinalsoc.eps','epsc')

% average computational time
figure('pos',[500 200 400 300])
boxplot(avgcomptime')
ylim([0 max(max(avgcomptime))*1.1])
% xlabel('case','FontSize',12)
ylabel('seconds','FontSize',12)
xticklabels({'unscheduled';'scheduled'})
saveas(gcf,'figs/avgcomptime.eps','epsc')


%% tight subplot

figure('pos',[500 200 400 300])
[ha, pos] = tight_subplot(1,3,[0 .05],[.1 .1],[.01 .01]) 
axes(ha(1));
boxplot(avgwait')
ylim([0 max(max(avgwait))*1.1])
ylabel('minutes','FontSize',12)
xticklabels({'U';'S'})
% title('average')
% text(0,0,'a')
title('a')

axes(ha(2));
boxplot(peakwait')
ylim([0 max(max(peakwait))*1.1])
ylabel('minutes','FontSize',12)
% xticklabels({'unscheduled';'scheduled'})
xticklabels({'U';'S'})
% title('peak')
% text(0,0,'b')
title('b')

axes(ha(3));
boxplot(totalcostfinalsoc'/1000)
hold on
weightedcosts=totalcostfinalsoc*elepprob;
plot(1:2,weightedcosts/1000,'xk')
offset1=0.2;
text(1+offset1,weightedcosts(1)/1000+offset1,sprintf('%0.0f yen',weightedcosts(1)))
text(2-3*offset1,weightedcosts(2)/1000+offset1,sprintf('%0.0f yen',weightedcosts(2)))
% ylim([0 max(max(totalcostfinalsoc))/1000*1.1])
ylabel('thousand yen','FontSize',12)
% xticklabels({'unscheduled';'scheduled'})
xticklabels({'U';'S'})
title('c')

% title('charging costs')

% eleptype='R';
% eleplist=[7416];
% eleplist=randi(10000,10,1);






% %% suplots
% 
% 
% figure('pos',[500 200 400 300])
% 
% subplot(1,3,1)
% boxplot(avgwait'/60)
% ylim([0 max(max(avgwait))/60*1.1])
% ylabel('minutes','FontSize',12)
% xticklabels({'unscheduled';'scheduled'})
% 
% subplot(1,3,2)
% boxplot(peakwait'/60)
% ylim([0 max(max(peakwait))/60*1.1])
% ylabel('minutes','FontSize',12)
% xticklabels({'unscheduled';'scheduled'})
% 
% subplot(1,3,3)
% boxplot(totalcostfinalsoc'/1000)
% hold on
% weightedcosts=totalcostfinalsoc*elepprob';
% plot(1:2,weightedcosts/1000,'xk')
% offset1=0.2;
% text(1+offset1,weightedcosts(1)/1000+offset1,sprintf('%0.0f yen',weightedcosts(1)))
% text(2-3*offset1,weightedcosts(2)/1000+offset1,sprintf('%0.0f yen',weightedcosts(2)))
% % ylim([0 max(max(totalcostfinalsoc))/1000*1.1])
% ylabel('thousand yen','FontSize',12)
% xticklabels({'unscheduled';'scheduled'})
% 
% 
% 
% %% suplots
% 
% 
% figure('pos',[500 200 400 300])
% 
% subplot(1,2,1)
% boxplot(avgwait'/60)
% hold on
% boxplot(peakwait'/60)
% ylim([0 max(max(peakwait))/60*1.1])
% ylabel('minutes','FontSize',12)
% xticklabels({'unscheduled';'scheduled'})
% 
% subplot(1,2,2)
% boxplot(totalcostfinalsoc'/1000)
% hold on
% weightedcosts=totalcostfinalsoc*elepprob';
% plot(1:2,weightedcosts/1000,'xk')
% offset1=0.2;
% text(1+offset1,weightedcosts(1)/1000+offset1,sprintf('%0.0f yen',weightedcosts(1)))
% text(2-3*offset1,weightedcosts(2)/1000+offset1,sprintf('%0.0f yen',weightedcosts(2)))
% % ylim([0 max(max(totalcostfinalsoc))/1000*1.1])
% ylabel('thousand yen','FontSize',12)
% xticklabels({'unscheduled';'scheduled'})
% 
% 
% %% suplots vertical
% 
% 
% figure('pos',[500 200 400 300])
% subplot(3,1,1)
% boxplot(avgwait'/60)
% ylim([0 max(max(avgwait))/60*1.1])
% ylabel('minutes','FontSize',12)
% set(gca,'xtick',[])
% 
% subplot(3,1,2)
% boxplot(peakwait'/60)
% ylim([0 max(max(peakwait))/60*1.1])
% ylabel('minutes','FontSize',12)
% set(gca,'xtick',[])
% 
% subplot(3,1,3)
% boxplot(totalcostfinalsoc'/1000)
% hold on
% weightedcosts=totalcostfinalsoc*elepprob';
% plot(1:2,weightedcosts/1000,'xk')
% offset1=0.2;
% text(1+offset1,weightedcosts(1)/1000+offset1,sprintf('%0.0f yen',weightedcosts(1)))
% text(2-3*offset1,weightedcosts(2)/1000+offset1,sprintf('%0.0f yen',weightedcosts(2)))
% % ylim([0 max(max(totalcostfinalsoc))/1000*1.1])
% ylabel('thousand yen','FontSize',12)
% xticklabels({'unscheduled';'scheduled'})




% %% tight subplot example
% 
% [ha, pos] = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01]) 
% for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end 
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
