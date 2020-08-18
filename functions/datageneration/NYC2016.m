
DataFolder=setDataFolder();
load([DataFolder 'scenarios/NYC2016'],'C')

for i=1:2
    FileIn=[DataFolder 'trips/NYC2016/d' num2str(i)];
    FileOut=[DataFolder 'trips/NYC2016-nodes/d' num2str(i)];
    StationRadius=[];
    WalkingSpeed=4;
    generateGPStrips(FileIn,FileOut,C,StationRadius,WalkingSpeed)
end
