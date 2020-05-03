%% 

function [T,C,A,Atimes,tripsubset]=generateScenario(Ac,Atimes,N,setlimits)


if nargin==4
    tripsubset=logical((Ac(:,1)>=setlimits(1,1)).*(Ac(:,1)<=setlimits(1,2)).*(Ac(:,2)>=setlimits(2,1)).*(Ac(:,2)<=setlimits(2,2)).* ...
                       (Ac(:,3)>=setlimits(1,1)).*(Ac(:,3)<=setlimits(1,2)).*(Ac(:,4)>=setlimits(2,1)).*(Ac(:,4)<=setlimits(2,2)));
else 
    tripsubset=true(length(Ac),1);
end

if 0
    figure
    scatter(A(tripsubset,1),A(tripsubset,2))
    axis equal
end


%% create C

Ac2=Ac(tripsubset,:);
[IDX,C,SUMD,D]=kmeans([Ac2(:,1:2);Ac2(:,3:4)],N);


%% create A, Atimes

A=reshape(IDX,sum(tripsubset),2);
Atimes=Atimes(tripsubset,:);


%% create T

origSelector=zeros(length(A),N);
destSelector=zeros(length(A),N);
for i=1:N
    origSelector(:,i)=logical(A(:,1)==i);
    destSelector(:,i)=logical(A(:,2)==i);
end
TravelTimes=min(60,Atimes(:,2)-Atimes(:,1));

T=zeros(N,N);
for i=1:N
    i
    for j=1:N
        Distances=TravelTimes(logical(origSelector(:,i).*destSelector(:,j)),:);
        if ~isempty(Distances)
            T(i,j)=mean(Distances);
        end
    end
end


%% reconstruct missing coordinates

addpath('../CarSharingModel/utilities/dijkstra_alg');
N=size(T,1);
for i=1:N
    for j=1:N
        if T(i,j)==0
            T(i,j)=dijkstra(T,i,j);
        end
    end
end


return



