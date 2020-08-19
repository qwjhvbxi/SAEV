
DataFolder=setDataFolder();
load([DataFolder 'scenarios/NYC2016'],'C')

for i=20:100
    i
    FileIn=[DataFolder 'trips/NYC2016/d' num2str(i)];
    FileOut=[DataFolder 'trips/NYC2016-nodes/d' num2str(i)];
    StationRadius=[];
    WalkingSpeed=4;
    generateGPStrips(FileIn,FileOut,C,StationRadius,WalkingSpeed)
end


%%

Ac=A;
Atimesc=Atimes;
Days=13:7:13+7*9;
for k=1:10
    A=Ac{k};
    Atimes=Atimesc{k};
    save([DataFolder 'trips/NYC2016/d' num2str(Days(k))],'A','Atimes');
end
