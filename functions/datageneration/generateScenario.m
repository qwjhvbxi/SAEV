%% [T,C,A,Atimes,tripsubset]=generateScenario(AcoordComplete,AtimesComplete,N[,setlimits,C])
% generate a new scenario and trip files from a larger scenario. The new
% scenario only considers the area within 'setlimits'
% [min_x,max_x,min_y,max_y] and with 'N' nodes 
% 
% generateScenario(AcoordComplete,AtimesComplete,N) 
%   generate new scenario with N nodes from same area 
% generateScenario(AcoordComplete,AtimesComplete,N,setlimits) 
%   generate new scenario with N nodes from area defined in setlimits
% generateScenario(AcoordComplete,AtimesComplete,N,setlimits,C) 
%   generate only trips from from area defined in setlimits and node
%   coordinates in C. N should be ==size(C,1)

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

    T=populateT(T);
    
else 
    
    T=[];
    
    Acoord=AcoordComplete(tripsubset,:);
    DistancesOrigin=(Acoord(:,1)*ones(1,N)-ones(length(Acoord),1)*C(:,1)').^2+(Acoord(:,2)*ones(1,N)-ones(length(Acoord),1)*C(:,2)').^2;
    DistancesDestination=(Acoord(:,3)*ones(1,N)-ones(length(Acoord),1)*C(:,1)').^2+(Acoord(:,4)*ones(1,N)-ones(length(Acoord),1)*C(:,2)').^2;
    [~,A1]=min(DistancesOrigin,[],2);
    [~,A2]=min(DistancesDestination,[],2);
    A=[A1 A2];
    
end

return



