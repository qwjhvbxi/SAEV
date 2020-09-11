
%% initializations

addpath functions utilities
DataFolder=setDataFolder();

Period=18:24;
gridfiles={'NY_DA_2016','Germany_DA_2019','Germany_DA_2019'};
gridoffsets=[0,-11,-11+7*20];    % 2016/01/01: Fri; 2019/01/01: Tue;

P=cpar('NYC2016');
P.Operations.maxwait=Inf; %20
P.e=1;
mvec=[10000 8000 6000 5000];

%% create Pmat

P.gridfile=gridfiles{1};
P.carbonprice=0;

P.enlayeralg='no';
P.m=mvec(1);
Pmat{1}=P;
P.m=mvec(2);
Pmat{2}=P;
P.m=mvec(3);
Pmat{3}=P;
P.m=mvec(4);
Pmat{4}=P;

P.enlayeralg='aggregate';
P.m=mvec(1);
Pmat{5}=P;
P.m=mvec(2);
Pmat{6}=P;
P.m=mvec(3);
Pmat{7}=P;
P.m=mvec(4);
Pmat{8}=P;

% P.gridfile=gridfiles{2};
% 
% P.gridfile=gridfiles{3};

GridOffsetsVec=repelem(gridoffsets,1,8);

%%

% load or create start variables to avoid saving files in parallel computing
load([DataFolder 'scenarios/' P.scenario],'T');
n=size(T,1);
StartFile=[DataFolder 'temp/StartFile_' num2str(n) '-' num2str(P.m) '.mat'];
if exist(StartFile,'file')
    load(StartFile,'StartSoc','StartPos');
else
    StartSoc=ones(1,P.m)*0.7;
    StartPos=randi(n,1,P.m);
    save(StartFile,'StartSoc','StartPos');
end

parfor k=1:length(Pmat)
    multiDaySim(Period,Pmat{k},GridOffsetsVec(k));
end

for k=1:length(Pmat)
    k
    [S1,~]=multiDaySim(Period,Pmat{k},GridOffsetsVec(k));
    S(k)=S1;
end


%% plots

plot([S(:).totavgwait])

plot(sum([S(:).cost]))
hold on


plot(prova,'.')



% k=2;
% % [S0,R]=multiDaySim(18,Pmat{k},GridOffsetsVec(k),1);
% [S1,~]=multiDaySim(Period,Pmat{k},GridOffsetsVec(k),0);
