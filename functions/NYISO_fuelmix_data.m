%% NYISO fuel mix data processing

addpath functions utilities
DataFolder=setDataFolder();

% Year='2018';
% MonthDays=[31 28 31 30 31 30 31 31 30 31 30 31];

Year='2016';
MonthDays=[31 29 31 30 31 30 31 31 30 31 30 31];

MonthDaysC=cumsum([0 MonthDays]);

%% load data

f1=[];

for m=1:12
    m
    for d=1:MonthDays(m)
        ThisDay=MonthDaysC(m)+d;
        Month=sprintf('%02d',m);
        Day=sprintf('%02d',d);
        FileName=[DataFolder 'input_files/NewYork/NYISO fuel mix/' Year Month Day 'rtfuelmix.csv'];
        g(ThisDay)=importdata(FileName);
        f1=[f1;g(ThisDay).data];
    end
end

Missing=0;
f2=[];
for k=1:length(g)
    
    L=length(g(k).data);
    Leff=min(L,Missing+24*12*7);
    f2=[f2;g(k).data(1:Leff)];
    Missing=max(0,k*24*12*7-length(f2));
end

f=f2;
Fuels=g(1).textdata(2:8,3);
GenerationByFuel_MW=reshape(f,length(Fuels),length(f)/length(Fuels))';

%% adjust change of order on 25 january 
% before 25 Jan:
% Dual Fuel
% Hydro
% Natural Gas
% Nuclear
% Other Fossil Fuels
% Other Renewables
% Wind
% ----
% from 25 Jan: 
% Dual Fuel
% Natural Gas
% Nuclear
% Other Fossil Fuels
% Other Renewables
% Wind
% Hydro
ChangeMoment=24*12*24; % 24 January 24:00
GenerationByFuel_MW_old=GenerationByFuel_MW;
GenerationByFuel_MW=[GenerationByFuel_MW_old(1:ChangeMoment,:);GenerationByFuel_MW_old(ChangeMoment+1:end,[1 7 2 3 4 5 6])];

% dual fuel mapped to natural gas, other renewables to biomass
CarbonEmissionsByFuel=[490, 490, 12, 700, 230, 11, 24]';
CarbonEmissions_TonPerMin=GenerationByFuel_MW*CarbonEmissionsByFuel/1000/60; % ton/minute
CarbonIntensity_kgPerMWh=CarbonEmissions_TonPerMin*1000*60./sum(GenerationByFuel_MW,2);

figure
plot(GenerationByFuel_MW)
figure
plot(CarbonEmissions_TonPerMin)
figure
plot(CarbonIntensity_kgPerMWh)

save([DataFolder 'grid/NYISO_fuelmix_' Year '-2.mat'],'Fuels','CarbonEmissionsByFuel','GenerationByFuel_MW','CarbonEmissions_TonPerMin','CarbonIntensity_kgPerMWh')


%% load data for prices

f1=[];

for m=1:12
    m
    for d=1:MonthDays(m)
        ThisDay=MonthDaysC(m)+d;
        Month=sprintf('%02d',m);
        Day=sprintf('%02d',d);
        FileName=[DataFolder 'input_files/NewYork/NYISO prices/' Year Month Day 'damlbmp_zone.csv'];
        g(ThisDay)=importdata(FileName);
        f1=[f1;g(ThisDay).data];
    end
end

f=f1(:,2);
Zones=g(1).textdata(2:16,2);
PriceByZone_USDMWh=reshape(f,length(Zones),length(f)/length(Zones))';

% NYC
% plot(PriceByZone_USDMWh(:,10))

save('data/eleprices/DayaheadPrices2016NY.mat','PriceByZone_USDMWh','Zones');


%% generate file with prices and carbon intensity
% results should be half-hourly

addpath functions utilities
DataFolder=setDataFolder();

load([DataFolder 'grid/NYISO_fuelmix_2016-2.mat'],'CarbonIntensity_kgPerMWh');
load('data/eleprices/DayaheadPrices2016NY.mat','PriceByZone_USDMWh');

DaysInYear=366;
HoursInYear=DaysInYear*24;

% prices
x=repelem(reshape(PriceByZone_USDMWh(1:HoursInYear),24,DaysInYear),2,1); % hourly to 2 per hour

% emissions
y1=mean(reshape(CarbonIntensity_kgPerMWh(1:HoursInYear*12),12/2,HoursInYear*2)); % 12 per hour to 2 per hour
y=reshape(y1',48,DaysInYear); 

save('data/eleprices/NY_DA_2016.mat','x','y');



