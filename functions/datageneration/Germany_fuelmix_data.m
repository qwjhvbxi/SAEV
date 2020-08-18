%% Germany fuel mix data processing

addpath functions utilities
DataFolder=setDataFolder();

% load([DataFolder 'eleprices/Generation2018Germany.mat'],'GermanyGeneration2018');
load([DataFolder 'eleprices/Generation2019Germany.mat'],'GermanyGeneration2019');
Generation=GermanyGeneration2019;

Fuels=Generation(:,2:end).Properties.VariableNames';
CarbonEmissionsByFuel=[230,820,490,820,650,0,0,24,24,12,700,230,45,230,11,11]'; % kg/MWh
GenerationByFuel_MW=table2array(Generation(:,2:end));                % MW (15 minutes)
CarbonEmissions_TonPerMin=GenerationByFuel_MW*CarbonEmissionsByFuel/1000/60;    % ton/minute
CarbonIntensity_kgPerMWh=CarbonEmissions_TonPerMin*1000*60./sum(GenerationByFuel_MW,2);

figure
plot(GenerationByFuel_MW)
figure
plot(CarbonEmissions_TonPerMin)
figure
plot(CarbonIntensity_kgPerMWh)

save([DataFolder 'grid/Germany_fuelmix_2019.mat'],'Fuels','CarbonEmissionsByFuel','GenerationByFuel_MW','CarbonEmissions_TonPerMin','CarbonIntensity_kgPerMWh')


%% generate file with prices and carbon intensity
% results should be half-hourly

addpath functions utilities
DataFolder=setDataFolder();

load([DataFolder 'grid/Germany_fuelmix_2019.mat'],'CarbonIntensity_kgPerMWh');
load([DataFolder 'eleprices/DayaheadPrices2018-2019Germany.mat'],'GermanyDayAheadEURMWh_2019');

% prices
x=repelem(reshape(GermanyDayAheadEURMWh_2019(1:8760),24,365),2,1); % hourly

% emissions
y1=mean(reshape(CarbonIntensity_kgPerMWh(1:(8760*4)),2,8760*2)); % quarter hourly
y=reshape(y1',48,365); 

save([DataFolder 'eleprices/Germany_DA_2019.mat'],'x','y');