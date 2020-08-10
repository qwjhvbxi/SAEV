
% % x1 x2 y1 y2
% % limits=[17 27 20 30];
% limits=[20 25 22 27];

function [OD,TPH,T]=generateTokyoScenario(scenarioname,limits,n,speedkmh,plots)


%% check for existing data

nomemap=['data/scenarios/' scenarioname '.mat'];
if exist(nomemap,'file')
    
    load(nomemap,'OD','TPH','T','gridcoord','pos','idx','clusters');
    
else
    
    
    %% load raw data
    
    load('data/sources/tripsCell.mat','Trips')
    load('data/sources/aev18data.mat','gridcoord','tph');

    
    %% resize map
    
    nostations=logical((gridcoord(:,1)<limits(1))+(gridcoord(:,2)<limits(3))+(gridcoord(:,1)>limits(2))+(gridcoord(:,2)>limits(4)));
    nostationsmapping=find(nostations);
    newstations=logical(1-nostations);
    newstationsmapping=find(newstations);
    oldstationsmapping=double(newstations);
    oldstationsmapping(newstations)=1:length(newstationsmapping);
    pos=gridcoord(newstations,:);
    for i=1:24
        Trips{i}=double(Trips{i});
        for j=1:length(nostationsmapping)
            Trips{i}(logical((Trips{i}(:,2)==nostationsmapping(j))+(Trips{i}(:,3)==nostationsmapping(j))),:)=[];
        end
        for j=1:length(newstationsmapping)
            Trips{i}(Trips{i}(:,2)==newstationsmapping(j),2)=oldstationsmapping(newstationsmapping(j));
            Trips{i}(Trips{i}(:,3)==newstationsmapping(j),3)=oldstationsmapping(newstationsmapping(j));
        end
        Trips{i}(:,1)=Trips{i}(:,1)/sum(Trips{i}(:,1));
    end
    
    
    %% clustering
    
    % find clusters
    [idx,clusters]=kmeans(pos,n);
    
    % find how many trips are between different zones
    tripsClustered=cell(24,1);              % all trips clustered 
    tripsClusteredInterzonal=cell(24,1);    % trips clustered between different zones
    interzonals=zeros(24,1); % fraction of trips between different zones
    for k=1:24 % for each hour
        tripsClustered{k}=[Trips{k}(:,1) idx(Trips{k}(:,2)) idx(Trips{k}(:,3))];
        interzonals(k)=(1-sum(tripsClustered{k}((tripsClustered{k}(:,2)-tripsClustered{k}(:,3)==0),1)));
        tripsClusteredInterzonal{k}=tripsClustered{k}( logical(1-(tripsClustered{k}(:,2)-tripsClustered{k}(:,3)==0))  , :);
        tripsClusteredInterzonal{k}(:,1)=tripsClusteredInterzonal{k}(:,1)/sum(tripsClusteredInterzonal{k}(:,1));
    end

    % calculate distance matrix
    delta=1.5; % tortuosity
    tempt=sqrt(  (clusters(:,1)*ones(1,n)-ones(n,1)*clusters(:,1)').^2 + (clusters(:,2)*ones(1,n)-ones(n,1)*clusters(:,2)').^2  );
    T=round(tempt*delta/speedkmh*60); % distance matrix (minutes)
    
    OD=tripsClusteredInterzonal;
    TPH=tph.*interzonals;
    TPH=TPH/sum(TPH);
    
    save(nomemap,'OD','TPH','T','gridcoord','pos','idx','clusters'); % last 4 are for plotting
    
end


if plots

    %% plot map of clusters
    figure
    scatter(pos(:,1),pos(:,2))
    hold on
    colori=hsv(10);
    for i=1:n
        scatter(pos((idx==i),1),pos((idx==i),2),[],ones(sum(idx==i),1)*i);
        scatter(clusters(i,1),clusters(i,2),[],i,'filled','o')
        text(clusters(i,1)+.3,clusters(i,2)+.1,int2str(i))
    end
    drawnow

    %% plot maps
    figure
    scatter(gridcoord(:,1),gridcoord(:,2))
    hold on
    scatter(pos(:,1),pos(:,2))
    axis equal

    figure
    scatter(pos(:,1),pos(:,2))
    axis equal

    pos(:,1)=pos(:,1)-limits(1);
    pos(:,2)=pos(:,2)-limits(3);

end

return

[OD,TPH,T]=generateTokyoScenario('TokyoSmall3',[20 25 22 27],10,30,0)







