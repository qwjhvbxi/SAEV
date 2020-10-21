

DataFolder=setDataFolder();
load([DataFolder 'scenarios/NYC2018'],'C')

for i=1:10
    FileIn=[DataFolder 'trips/NYC2018_10wed_small/d' num2str(i)];
    FileOut=[DataFolder 'trips/NYC2018_10wed_small_nodes/d' num2str(i)];
    StationRadius=[];
    WalkingSpeed=4;
    generateGPStrips(FileIn,FileOut,C,StationRadius,WalkingSpeed)
end

