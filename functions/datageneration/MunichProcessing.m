
DataFolder=setDataFolder();

% odagents
% uiopen([DataFolder 'input_files/Munich/od_agents.csv',1)

A=table2array(odagents(:,1:2))+1;
TimeVec=table2array(odagents(:,4));
Atimes=repmat(TimeVec.Hour*60+TimeVec.Minute,1,2)+1;

Adist=table2array(odagents(:,3));

n=max(A(:));

Tmeters=sparse(n,n);

%%

for i=1:n
    Aind=find(A(:,1)==i);
    An=A(Aind,2);
    for j=1:length(An)
        if Tmeters(i,An(j))==0
            Tmeters(i,An(j))=Adist(Aind(j));
        else
            if Tmeters(i,An(j))~=Adist(Aind(j))
                Tmeters(i,An(j))
                Adist(Aind(j))
            end
        end
    end
end

T=round(Tmeters/1000/30*60);

[A,Atimes,~]=cleanData(A,Atimes);

% remove trips inside same node
SameNode=logical(A(:,1)==A(:,2));
A(SameNode,:)=[];
Atimes(SameNode,:)=[];

save([DataFolder 'trips/Munich.mat'],'A','Atimes')
save([DataFolder 'scenarios/Munich.mat'],'T');



