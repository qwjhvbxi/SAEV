

%% individual results plots
% use with files in 'out' folder

convertvarformat



% total SOC
figure('pos',[500 200 400 300])
plot(linspace(0,24,P.tsim+1),q_check')
xlim([0 24])
% ylim([0.6 0.9])
xticks((0:4:24))
xlabel('hour','FontSize',12)
saveas(gcf,'figs/soc1.eps','epsc')


% price electricity and charge power
figure('pos',[500 200 400 300])%(100)
% plot(linspace(0,24,P.tsim),P.elep(1:P.tsim))
area(linspace(0,24,P.tsim),P.elep(1:P.tsim),'FaceColor',[.6,.8,1],'LineStyle','none')
ylabel('price [JPY/kWh]')
hold on
yyaxis right
plot(linspace(0,24,P.tsim),(sum(e_check)-sum(e2_check))*P.battery*60/2,'r') 
set(gca,'YColor','k');
legend({'price of electricity';'charge power'})
xlim([0 24])
xticks((0:4:24))
ylabel('charge power [kW]')
xlabel('hour')
saveas(gcf,'figs/soc1.eps','epsc')



return


%% other plots

% waiting people
figure
plot(linspace(0,24,P.tsim+1),squeeze(sum(d_check,2))')
xlim([0 24])
xticks((0:4:24))
xlabel('hour')

figure
plot(P.elep(1:P.tsim))
hold on
yyaxis right
plot(e_check','-')
legend({'price of electricity';'energy exchange'})






%% length of simulation
figure
plot(S.durata)
title('length of simulation micro')
figure
plot(S.duratamacro)
title('length of simulation macro')

% status total
figure
colormap(gray(3)) 
im=squeeze(sum(u_check(:,:,:))) + squeeze(sum(sum(p_check(:,:,:,:)),3))*2; % 1 = waiting, 2 = moving
image(im)


%% figures for specific car

car=1;

% charge and electricity prices
figure
plot(linspace(0,24,1440/P.e),e_check(car,:))
ylabel('power [kW]');
hold on
yyaxis right
plot(linspace(0,24,1440/P.e),P.elep(1:P.tsim))
ylabel('price [$/kWh]');
legend({'energy exchange';'price of electricity'})
axis([0 24 0 1])
xticks((0:4:24))
xlabel('hour')

% SOC
figure
plot(linspace(0,24,1440/P.e+1), q_check(car,:))
title('SOC')
axis([0 24 0 1])
xticks((0:4:24))
xlabel('hour')

% status
figure
subplot(3,1,1)
area(squeeze(sum(sum(p_check(:,car,:,:)))));
title('moving')
subplot(3,1,2)
area(squeeze(sum(u_check(:,car,:))))
title('at station')
subplot(3,1,3)
area(econ_check(car,:))
title('can''t move')


%% map and charging stations for Tokyo

% plot map of clusters
figure
% limits=[17 27 20 30];
limits=[20 25 22 27];
% topright=[35.824811, 139.957411];

% calculate coordinates of limits

% nometrips=['data/tripsClustered' P.mapmod '-' num2str(P.n)];
nometrips=['data/tripsClusteredSmall-10'];
load(nometrips)

topright=[35.824811, 139.957411];
lat1=topright(1)*pi/180;
lon1=topright(2)*pi/180;


x1=-(40-limits(2))/6371;
y1=-(40-limits(4))/6371;

lon2=lon1+x1/cos(lat1);
lat2=lat1+y1;
        
lat2=lat2/pi*180;
lon2=lon2/pi*180;

[lat2 ;lon2]

x2=-(40-limits(1))/6371;
y2=-(40-limits(3))/6371;

lon3=lon1+x2/cos(lat1);
lat3=lat1+y2;
        
lat3=lat3/pi*180;
lon3=lon3/pi*180;

[lat3 ;lon3]

DataFolder=setDataFolder();
% M=imread('data/TokyoStreetMapNarrow.png');
% M=imread('data/TokyoStreetMapSmall2.png');
M=imread([DataFolder 'input_files\Tokyo\maps/TokyoStreetMap5km.png']);
% image([0 10],[0 10],flipud(M)/1.5)
% image([0 10],[0 10],254/2+flipud(M)/2)
image([0 5],[0 5],254/2+flipud(M)/2)
set(gca,'ydir','normal');

hold on

% k=P.n;
k=10;
% colori=hsv(k);
colori=colorcube(k*2)/1.5;
% colori(sum(colori,2)>0.8,:)=[];
% colori=lines(k);
for i=1:k
    
    % small node
    scatter(gridcoord((idx==i),1),gridcoord((idx==i),2),20,ones(sum(idx==i),1)*colori(i,:),'filled','o');

    % cluster centroids
    scatter(clusters(i,1),clusters(i,2),150,colori(i,:),'filled','o')
    
    % centroid numbers
    text(clusters(i,1)+.3,clusters(i,2)+.1,int2str(i),'Color','k')
%     text(clusters(i,1),clusters(i,2),int2str(i),'Color','w')

end

xlabel('km')
ylabel('km')
% scatter(Param.gridcoord(Param.charginst,1),Param.gridcoord(Param.charginst,2),'r','filled');
% scatter(Param.charginst(:,1),Param.charginst(:,2),'r','filled');
axis([0 limits(2)-limits(1) 0 limits(4)-limits(3)])
axis square;
