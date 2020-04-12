
%% compare results with aggregate
% use with files in 'out' folder

convertvarformat

% power exchanged theory/effective
ptheo=((charge_check(1,:)-charge_check(2,:)).*av_check)'/P.macrosteps/P.e*60; % kWh, av_check kWh!!!
peffe=sum(e_check-e2_check)'*P.battery/P.e*60; 
peffeavg=repelem(sum(reshape(peffe,P.macrosteps,P.tsim/P.macrosteps))/P.macrosteps,1,P.macrosteps)';

% SOC theory/effective
stheo=(cumsum(ptheo)-etrip_check')/P.battery/P.m+P.initialsoc;
seffe=mean(q_check(:,2:end))';

% lower layer stats
moving=squeeze(sum(sum(p_check(:,:,:,:)),3));
notmoving=sum(1-moving)';
dktripeff=sum(reshape(sum(moving(:,2:end)),[P.macrosteps,P.T]))';


%% plots

% energy availability / theory vs practice
figure('pos',[500 200 400 300])
plot(linspace(0,24,P.T),max(0,P.m*P.macrosteps-[P.dktrip(1:P.T) dktripeff]))
xlim([0,24])
xticks((0:4:24))
xlabel('hour','FontSize',12)
ylabel('connected vehicles-steps','FontSize',12)
legend({'energy layer';'transport layer'},'FontSize',12)
saveas(gcf,'figs/connected.eps','epsc')

% total SOC
figure('pos',[500 200 400 300])
plot(linspace(0,24,P.tsim+1),q_check')
xlim([0 24])
xticks((0:4:24))
xlabel('hour','FontSize',12)
saveas(gcf,'figs/soc1.eps','epsc')

% total SOC as image
figure('pos',[500 200 400 300])
colormap(flipud(jet))
imagesc([0 24],[0 P.m],q_check)
caxis([0 1])
colorbar
xlim([0 24])
xticks((0:4:24))
xlabel('hour','FontSize',12)
ylabel('vehicle','FontSize',12)
saveas(gcf,'figs/soc2.eps','epsc')

% status as image
% status: 1 = waiting, 2 = moving, 3 = charging, 4 = v2g
figure('pos',[500 200 400 300])
themap=[0.2 0.1 0.5
        %0.1 0.5 0.8
        0.2 0.7 0.6
        0.8 0.7 0.3
        0.9 1 0];
colormap(themap) 
im=     squeeze(sum(u_check(:,:,:))) + ...
        squeeze(sum(sum(p_check(:,:,:,:)),3))*2 + ...
        [zeros(P.m,1) (e_check>0)*2] + ...
        [zeros(P.m,1) (e2_check>0)*3] ;
image([0 24],[0 P.m],im)
cb=colorbar;
cb.Ticks=[1.5:5.5];
cb.TickLabels={'idle','moving','charging','v2g'};
xlabel('hour','FontSize',12)
ylabel('vehicle','FontSize',12)
saveas(gcf,'figs/statuses.eps','epsc')

% power exchanged
figure('pos',[500 200 400 300])
area(linspace(0,24,P.tsim),P.elep(1:P.tsim),'FaceColor',[190,220,210]/255,'LineStyle','none')
ylabel('price [JPY/kWh]','FontSize',12)
hold on
yyaxis right
plot(linspace(0,24,P.tsim),ptheo/1000,'r-')
plot(linspace(0,24,P.tsim),peffe/1000,'k-')
% plot(linspace(0,24,P.tsim),peffeavg/1000,'r--')
set(gca,'YColor','k');
legend({'el. price';'energy layer';'transp. layer'},'Location','NorthOutside','Orientation','horizontal','FontSize',10)
xlim([0 24])
xticks((0:4:24))
ylabel('charge power [MW]','FontSize',12)
xlabel('hour','FontSize',12)
saveas(gcf,'figs/power.eps','epsc')



if 0

% soc average
figure
plot([stheo seffe])

% available vs power exchange
figure
area([notmoving(2:end)*P.ac*P.battery ],'LineStyle','None')
hold on
plot(abs(peffe),'k-','LineWidth',1.5)

prod(notmoving(2:end)*P.ac*P.battery >= abs(peffe))==1 % check if correct (1)
sum(sum((e_check>0).*(e2_check>0))) %check if charging and discharging at same time (0)
sum(sum((charge_check(1,:)>0).*(charge_check(2,:)>0)))

end



