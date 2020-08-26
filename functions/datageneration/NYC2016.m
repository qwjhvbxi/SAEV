
DataFolder=setDataFolder();
load([DataFolder 'scenarios/NYC2016'],'T','C');

for i=39:80
    i
    FileIn=[DataFolder 'trips/NYC2016/d' num2str(i)];
    FileOut=[DataFolder 'trips/NYC2016-nodes/d' num2str(i)];
    StationRadius=[];
    WalkingSpeed=4;
    generateGPStrips(FileIn,FileOut,C,StationRadius,WalkingSpeed)
end

%% generate trip stats

DataFolder=setDataFolder();
load([DataFolder 'scenarios/NYC2016'],'T','C');

for i=37:60
    i
    TripFile=[DataFolder 'trips/NYC2016-nodes/d' num2str(i)];
    TripName='NYC2016-nodes';
    EMDFileName=[TripName '-' num2str(i)];
    etsim=48;
    load(TripFile,'A','Atimes');
    generateEMD(A,Atimes,T,etsim,EMDFileName)
end



%%
% 
% Ac=A;
% Atimesc=Atimes;
% Days=13:7:13+7*9;
% for k=1:10
%     A=Ac{k};
%     Atimes=Atimesc{k};
%     save([DataFolder 'trips/NYC2016/d' num2str(Days(k))],'A','Atimes');
% end
