function DataFolder=setDataFolder()

if strcmp(version('-release'),'2019b')
    DataFolder='c:\Users\Riccardo\Documents\model_aev5_data\';
    return
end
if strcmp(version('-release'),'2015a')
    DataFolder='data/';
    return
end
DataFolder='c:\Users\iacob\Documents\model_aev5_data\';

% if on cluster, change