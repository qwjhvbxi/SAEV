



%% plot results

scenarioname='DayaheadPrices2017Germany';   % germany
% scenarioname='DayaheadPrices2016Tokyo';     % tokyo
load([DataFolder , scenarioname '-r.mat'])
resn=[21:104];              % which files to load (in case all)


resn=[104:144];              % which files to load (in case all)


z=outdesc(resn);            % load files


if 0        % plot profiles selected

    % germany
    load([DataFolder , scenarioname])
    phour=DayaheadPrices2017Germany;
    x=reshape(phour,24,365);

    % tokyo
    load('data/JEPX');
    TokyoDAFY2016=(Jepx.M((365*3+1)*48+1:end,8));
    x=reshape(TokyoDAFY2016,48,365);

    figure
    stairs(0:0.5:23.5,x(:,dayssel));
    % legend(cellstr(num2str(dayssel)))
    xlim([0,24])
    xlabel('hour')
    ylabel('yen/kWh')

end


% analyze results
% dayssel=[u(1:10)];    % select which eleprofiles to plot
% dayssel=[u([1:7,9:10])];    % select which eleprofiles to plot
seltypes='R';             % which types to plot ... not used now






y=z(1);     % initialize
for k=1:length(seltypes)
    for i=1:length(dayssel)
        for j=1:4

            y(i,j)=z(find(([z.type]==j).*([z.eleprofseed]==dayssel(i)) ...
                ....*([z.eleproftype]=='T')...
                ,1));

        end
    end
end

% create matrices
w=zeros(9,length(dayssel),4);  
for i=1:4
    tc(i,:)=([y(1:end,i).totalcostfinalsoc]); % total costs
    w(:,:,i)=[y(1:end,i).waittimes]; % waiting times
    comp(:,:,i)=[y(1:end,i).comptime]; % computational times
end

mean(comp,2)


% plot waiting times
figure('pos',[10 10 400 400])
hold on
for i=1:4
plot(squeeze(w(9,:,i))/60,markervec(i,1))
end
xticklabels(dayssel);
% title('peak wait times (min)')
ylabel('peak wait times (min)')
ylabel('case')
legend({'unscheduled';'1-layer';'2-layer';'2-layer +V2G'})


%% boxplots


% avgpeak waiting  

figure('pos',[1000 200 250 250])
hold on
boxplot(squeeze(w(7,:,:))/60)
ylabel('peak wait time (min)','FontSize',12)
xlabel('case','FontSize',12)
% ylim([0,7])
% xticklabels({'unscheduled';'1-layer';'2-layer';'2-layer +V2G'})
saveas(gcf,'figs/boxplotwaitpeakrealR.eps','epsc')

% peak waiting  

figure('pos',[1000 200 250 250])
hold on
boxplot(squeeze(w(9,:,:))/60)
ylabel('peak wait times (min)','FontSize',12)
xlabel('case','FontSize',12)
% ylim([0,7])
% xticklabels({'unscheduled';'1-layer';'2-layer';'2-layer +V2G'})
saveas(gcf,'figs/boxplotwaitpeakR.eps','epsc')

% costs 

figure('pos',[1000 200 250 250])
% boxplot(tc'/mean([y.totalkm]))
% ylabel('cost/km (yen)','FontSize',12)
boxplot(tc')
ylabel('total cost (yen)','FontSize',12)
xlabel('case','FontSize',12)
saveas(gcf,['figs/boxplottotalcost' seltypes '.eps'],'epsc')


%% plot electricity prices

seltypes='R';             % which types to plot ... not used now
dayssel= [6713
        7153
        6421
        4191
        3908
        8162
        3175
        8146
        7891
        8523];

% seltypes='T';
% dayssel=[356,127,8,188,284,185,3,26,274,138];

prices2=zeros(48,length(dayssel));
for i=1:length(dayssel)
    
    load(['../data/elep/elep-' seltypes  num2str(dayssel(i))])
    prices2(:,i)=elep(1:48);

end

%% stairs

figure('pos',[500 200 400 300])
stairs(0:0.5:23.5,prices);
xlim([0,24])
xlabel('hour')
ylabel('yen/kWh')
box on


%% heat maps 

figure('pos',[500 200 400 300])
colormap jet(30)
image(0:0.5:23.5,1:10,prices')
yticks(1:10)
xlabel('hour')
ylabel('profile')
h = colorbar;
ylabel(h, 'yen/kWh')

figure('pos',[500 200 400 300])
colormap jet(100)
image(0:0.5:23.5,1:10,prices')
yticks(1:10)
xlabel('hour')
ylabel('profile')
h = colorbar;
ylabel(h, 'yen/kWh')

figure('pos',[500 200 400 300])
colormap jet(100)
image(0:0.5:23.5,1:10,prices2')
yticks(1:10)
xlabel('hour')
ylabel('profile')
h = colorbar;
ylabel(h, 'yen/kWh')
% box on


%% stairs vecchio metodo

figure('pos',[500 200 400 300])
hold on
for i=1:length(dayssel)
    
    load(['../data/elep/elep-' seltypes  num2str(dayssel(i))])
    stairs(0:0.5:23.5,elep(1:48));

end
xlim([0,24])
xlabel('hour')
ylabel('yen/kWh')
box on



% random profile selected ID: 6421
outdesc([107 121:123])


figure
hold on
plot(tc,'o')

% normalized results
normalizationvec=ones(10,1)*0.1;
normalizationvec=r(1:10,10)
% normalizationvec=r([1:7,9:10],10);
tc_norm=tc*normalizationvec
pw_norm1=squeeze(w(7,:,:))'/60*normalizationvec
pw_norm2=squeeze(w(9,:,:))'/60*normalizationvec

figure
plot(tc_norm)
hold on
yyaxis right
plot(pw_norm)






figure
hold on
for i=1:4
plot(squeeze(w(9,:,i))/60,markervec(i,1))
end
xticklabels(dayssel);
% title('peak wait times (min)')
ylabel('peak wait times (min)')
ylabel('case')
legend({'unscheduled';'1-layer';'2-layer';'2-layer +V2G'})






 
figure
plot(wt')

figure
plot(wtavg')
xticklabels(eleprofs);
title('average wait times (sec)')

figure
plot(tc')
xticklabels(eleprofs);
title('total cost (yen) normalized with final soc')







 
% q=reducedscenarios(:,2)
% [a,b]=sort(q,'descend');
% [u(b) q(b) cumsum(q(b))]
% 
% 
% % elep profiles
% elepday=zeros(n,1);
% for i=1:n
%     elepday(i)=str2double(z(i).eleprofseed(2:end));
% end
% % eleprofs=unique(elepday);
% eleprofs=elepday(1:4:end);


% waittimes=[z.waittimes];
% 
% wt=reshape(waittimes(9,:),4,n/4);
% 
% wtavg=reshape(waittimes(8,:),4,n/4);
% 
% reshape([z.type],4,n/4)
% 
% 
% eleprofsp=zeros(length(eleprofs),5);
% for i=1:length(eleprofs)
%     eleprofsp(i,:)=reducedscenarios((reducedscenarios(:,1)==eleprofs(i)),:);
% end
% 
% 
% ct0=reshape([z.totalcost],4,n/4);
% ct=reshape([z.totalcostfinalsoc],4,n/4);
% et=reshape(elepday,4,n/4);
% 
% 
% 
% figure
% plot(wt'/60)
% xticklabels(eleprofs);
% title('peak wait times (min)')
% 
% figure
% plot(wtavg')
% xticklabels(eleprofs);
% title('average wait times (sec)')
% 
% figure
% plot(ct0')
% xticklabels(eleprofs);
% title('total cost (yen)')
% 
% figure
% plot(ct')
% xticklabels(eleprofs);
% title('total cost (yen) normalized with final soc')
% 
% 
% wt*eleprofsp(:,4)
% wtavg*eleprofsp(:,4)
% ct0*eleprofsp(:,4)
% 
% plot(ct*eleprofsp(:,4))
% ylim([0,4*10^4])
% 
% 
% figure
% hold on
% for i=1:9
% plot(1:4,ct(:,i),'LineWidth',log(eleprofsp(i,4)+1)*1000)
% % xticklabels(eleprofs);
% end
% title('total cost (yen) normalized with final soc')















