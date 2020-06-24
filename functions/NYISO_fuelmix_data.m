% NYISO fuel mix data processing
DataFolder=setDataFolder();


Year='2018';
MonthDays=[31 28 31 30 31 30 31 31 30 31 30 31];

% Year='2016';
% MonthDays=[31 29 31 30 31 30 31 31 30 31 30 31];

MonthDaysC=cumsum([0 MonthDays]);
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

% dual fuel mapped to natural gas, other renewables to biomass
CarbonEmissionsByFuel=[490, 490, 12, 700, 230, 11, 24]';
Fuels=g(1).textdata(2:8,3);
GenerationByFuel_MW=reshape(f,length(Fuels),length(f)/length(Fuels))';
CarbonEmissions_TonPerMin=GenerationByFuel_MW*CarbonEmissionsByFuel/1000/60; % ton/minute
CarbonIntensity_kgPerMWh=CarbonEmissions_TonPerMin*1000*60./sum(GenerationByFuel_MW,2);

figure
plot(GenerationByFuel_MW)
figure
plot(CarbonEmissions_TonPerMin)
figure
plot(CarbonIntensity_kgPerMWh)

save([DataFolder 'grid/NYISO_fuelmix_' Year '.mat'],'Fuels','CarbonEmissionsByFuel','GenerationByFuel_MW','CarbonEmissions_TonPerMin','CarbonIntensity_kgPerMWh')

