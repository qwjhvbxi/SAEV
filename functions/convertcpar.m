function convertcpar(Scenario)
DataFolder=getdatafolder();
P=cpar(Scenario);
save([DataFolder 'par/' Scenario '.mat'],'P');
