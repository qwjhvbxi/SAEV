
addpath functions utilities
DataFolder=setDataFolder();


%% non-clustered

if 0
    
    % uiopen([DataFolder 'input_files/Munich/od_agents.csv',1)

    A=table2array(odagents(:,1:2))+1;
    TimeVec=table2array(odagents(:,4));
    Adist=table2array(odagents(:,3));

end


%% clustered

% open file

i=1;
% FileID='';
FileID=['_' num2str(i)];
% odagentsclustered=readtable([DataFolder 'input_files/Munich/od_agents_clustered' FileID '.csv'],'DatetimeType','datetime');

%% 

Atot=[];
Aduratot=[];
Adisttot=[];

for i=1:7

    % terrible but fast solution
    eval(['odac=odagentsclustered' num2str(i) ';'])
    % odac=odagentsclustered1;

    A=table2array(odac(:,1:2))+1;
    TimeVec=table2array(odac(:,3));
    Adura=table2array(odac(:,4));
    Adist=table2array(odac(:,5));
    
    Atot=[Atot;A];
    Aduratot=[Aduratot;Adura];
    Adisttot=[Adisttot;Adist];

    Atimes=repmat(TimeVec.Hour*60+TimeVec.Minute,1,2)+1;

    [A,Atimes,~]=cleanData(A,Atimes,1);
    save([DataFolder 'trips/Munich_1week/d' num2str(i) '.mat'],'A','Atimes')

end

%% calculations




%% create matrix T

n=max(Atot(:));
Tmeters=zeros(n,n);
Tseconds=zeros(n,n);
N=zeros(n,n);

for i=1:n
    Aind=find(Atot(:,1)==i);
    An=Atot(Aind,2);
    for j=1:length(An)
        Tmeters(i,An(j))=Tmeters(i,An(j))+Adisttot(Aind(j));
        Tseconds(i,An(j))=Tseconds(i,An(j))+Aduratot(Aind(j));
        N(i,An(j))=N(i,An(j))+1;
%         if Tmeters(i,An(j))==0
%             Tmeters(i,An(j))=Adist(Aind(j));
%             Tseconds(i,An(j))=Adura(Aind(j));
%         else
%             if Tmeters(i,An(j))~=Adist(Aind(j))
%                 Tmeters(i,An(j))
%                 Adist(Aind(j))
%                 Adura(Aind(j))
%             end
%         end
    end
end

% T=round(Tmeters/1000/30*60);
T=round(Tseconds./N/60);
T(isnan(T))=0;

T=populateT(T);

T(1:n+1:end)=0;



% save([DataFolder 'trips/Munich.mat'],'A','Atimes')
% save([DataFolder 'scenarios/Munich.mat'],'T');
save([DataFolder 'scenarios/Munich_clustered_week.mat'],'T');



