%% 

function [T,C,A,Atimes,tripsubset]=generateScenario(AcoordComplete,AtimesComplete,N,setlimits,C)

if nargin>=4 && ~isempty(setlimits)
    tripsubset=logical((AcoordComplete(:,1)>=setlimits(1,1)).*(AcoordComplete(:,1)<=setlimits(1,2)).* ...
                       (AcoordComplete(:,2)>=setlimits(2,1)).*(AcoordComplete(:,2)<=setlimits(2,2)).* ...
                       (AcoordComplete(:,3)>=setlimits(1,1)).*(AcoordComplete(:,3)<=setlimits(1,2)).* ...
                       (AcoordComplete(:,4)>=setlimits(2,1)).*(AcoordComplete(:,4)<=setlimits(2,2)));
else 
    tripsubset=true(length(AcoordComplete),1);
end

if 0
    figure
    scatter(A(tripsubset,1),A(tripsubset,2))
    axis equal
end

%% create Atimes

Atimes=AtimesComplete(tripsubset,:);


if nargin<5

    %% create C

    Acoords=AcoordComplete(tripsubset,:);
    [IDX,C,SUMD,D]=kmeans([Acoords(:,1:2);Acoords(:,3:4)],N);


    %% create A

    A=reshape(IDX,sum(tripsubset),2);


    %% create T

    TravelTimes=min(60,Atimes(:,2)-Atimes(:,1));

    % fix travel times for trips ending next day
    TravelTimes(TravelTimes<-1380)=TravelTimes(TravelTimes<-1380)+1440;

    % remove errors/outliers
    ProblemTrips=(TravelTimes<0);
    TravelTimes(ProblemTrips)=[];
    A(ProblemTrips,:)=[];
    Atimes(ProblemTrips,:)=[];

    origSelector=zeros(length(A),N);
    destSelector=zeros(length(A),N);
    for i=1:N
        origSelector(:,i)=logical(A(:,1)==i);
        destSelector(:,i)=logical(A(:,2)==i);
    end


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

else 
    
    T=[];
    
    Acoord=AcoordComplete(tripsubset,:);
    Distances=(Acoord(:,1)*ones(1,N)-ones(length(Acoord),1)*C(:,1)')^2+(Acoord(:,2)*ones(1,N)-ones(length(Acoord),1)*C(:,2)')^2;
    [~,A]=min(Distances,[],2);
    
end

return



