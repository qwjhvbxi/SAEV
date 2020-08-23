
addpath functions
DataFolder=setDataFolder();

%% clustered

% uiopen([DataFolder 'input_files/Munich/od_agents_clustered.csv',1)

A=table2array(odagentsclustered(:,1:2))+1;
TimeVec=table2array(odagentsclustered(:,3));
Adura=table2array(odagentsclustered(:,4));
Adist=table2array(odagentsclustered(:,5));



if 0

    %% odagents
    
    % uiopen([DataFolder 'input_files/Munich/od_agents.csv',1)

    A=table2array(odagents(:,1:2))+1;
    TimeVec=table2array(odagents(:,4));
    Adist=table2array(odagents(:,3));

end

%% calculations

Atimes=repmat(TimeVec.Hour*60+TimeVec.Minute,1,2)+1;
n=max(A(:));
Tmeters=zeros(n,n);
Tseconds=zeros(n,n);
N=zeros(n,n);

for i=1:n
    Aind=find(A(:,1)==i);
    An=A(Aind,2);
    for j=1:length(An)
        Tmeters(i,An(j))=Tmeters(i,An(j))+Adist(Aind(j));
        Tseconds(i,An(j))=Tseconds(i,An(j))+Adura(Aind(j));
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
T(1:n+1:end)=0;
T(isnan(T))=0;

T=populateT(T);

[A,Atimes,~]=cleanData(A,Atimes,1);



% save([DataFolder 'trips/Munich.mat'],'A','Atimes')
% save([DataFolder 'scenarios/Munich.mat'],'T');
save([DataFolder 'trips/Munich_clustered.mat'],'A','Atimes')
save([DataFolder 'scenarios/Munich_clustered.mat'],'T');



