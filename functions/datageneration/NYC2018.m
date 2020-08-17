

DataFolder=setDataFolder();
load('data/scenarios/NYC2018','C')

for i=1:2
    FileIn=[DataFolder 'trips/NYC2018_10wed/d' num2str(i)];
    FileOut=[DataFolder 'trips/NYC2018_10wed-nodes/d' num2str(i)];
    StationRadius=[];
    WalkingSpeed=4;
    generateGPStrips(FileIn,FileOut,C,StationRadius,WalkingSpeed)
end

