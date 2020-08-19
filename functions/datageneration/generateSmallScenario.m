%% generate smaller scenarios
% generate C coordinates: one day is ok, with optional limits (this file)
% generate T times: need as much data as possible
% generate tripfile: need C and T, optional limits

%% generate scenario nodes and travel times files

if 0
    
    addpath functions
    DataFolder=setDataFolder();
    P=cpar('NYC2016');

    setlimits=[0,5;0,5];
    
    tripFileLocation=[DataFolder 'trips/' P.tripfile '.mat'];
    load(tripFileLocation,'A','Atimes');
    Ac=A;
    Actimes=Atimes;

    [T,C,A2,A2times,tripsubset]=generateScenario(Ac{P.tripday},Actimes{P.tripday},10,setlimits);

    % plot trips coordinates
    Ac10=Ac{P.tripday}(tripsubset,:);
    colori=lines(10);
    figure
    hold on
    for i=1:10
        scatter(Ac10(A2(:,1)==i,1),Ac10(A2(:,1)==i,2),1,colori(i,:))
    end

    % plot nodes coordinates
    scatter(C(:,1),C(:,2),5)
    axis equal

    % save C and T
    save([DataFolder 'scenarios/NYC2016-small.mat'],'C','T');

end


%% generate trips file

addpath functions
DataFolder=setDataFolder();
load([DataFolder 'scenarios/NYC2016-small.mat'],'C','T');

setlimits=[0,5;0,5];

if 0  % specific day

    P=cpar('NYC2016');
    tripFileLocation=[DataFolder 'trips/' P.tripfile '.mat'];

    % load original files
    load(tripFileLocation,'A','Atimes');
    Ac=A;
    Actimes=Atimes;
    Ac10=Ac{P.tripday}(tripsubset,:);

    [~,~,A2,A2times,tripsubset]=generateScenario(Ac{P.tripday},Actimes{P.tripday},10,setlimits,C);

    % remove trips with same origin/destination
    DifferentOD=(A2(:,1)~=A2(:,2));
    A3=A2(DifferentOD,:);
    A3times=A2times(DifferentOD,:);
    Ac10_2=Ac10(DifferentOD,:);

    k=300;
    [A,Atimes]=reduceTrips(A3,A3times,k)

    save([DataFolder 'trips/NYC2016-small_13Jan.mat'],'A','Atimes');

end



if 1  % folder based
    
    TripRatio=1/100;%1/200; % 
    Period=1:76;
    TripFolder='NYC2016';
    NewTripFolder='NYC2016-small_100';
    
    for d=1:length(Period)
        
        tripFileLocation=[DataFolder 'trips/' TripFolder '/d' num2str(Period(d)) '.mat'];
        load(tripFileLocation,'A','Atimes');

        [~,~,A2,A2times,tripsubset]=generateScenario(A,Atimes,10,setlimits,C);
        
        [A3,A3times,~]=cleanData(A2,A2times,1);

%         % remove trips with same origin/destination
%         DifferentOD=(A2(:,1)~=A2(:,2));
%         A3=A2(DifferentOD,:);
%         A3times=A2times(DifferentOD,:);

        k=round(length(A3)*TripRatio);
        [A,Atimes]=reduceTrips(A3,A3times,k);

        save([DataFolder 'trips/' NewTripFolder '/d' num2str(Period(d)) '.mat'],'A','Atimes');
    end
    
end



%% even smaller scenario
if 0

    load([DataFolder 'trips/NYC2016-small_13Jan.mat'],'A','Atimes');

    k=100;
    [A,Atimes]=reduceTrips(A,Atimes,k)
    save([DataFolder 'trips/NYC2016-small2_13Jan.mat'],'A','Atimes');

end




function [A,Atimes]=reduceTrips(A3,A3times,k)

M=sort(randperm(length(A3),k));
A=A3(M,:);
Atimes=A3times(M,:);

% plot
if 0
    Ac10_3=Ac10_2(M,:);
    colori=lines(10);
    figure
    hold on
    for i=1:10
        scatter(Ac10_3(A(:,1)==i,1),Ac10_3(A(:,1)==i,2),5,colori(i,:))
    end
    scatter(C(:,1),C(:,2),50)
    axis equal

    figure
    plot(histc(Atimes(:,1),0:10:1440))

end

end