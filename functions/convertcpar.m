function convertcpar(Scenario)
DataFolder=setDataFolder();
P=cpar(Scenario);
save([DataFolder 'par/' Scenario '.mat'],'P');
