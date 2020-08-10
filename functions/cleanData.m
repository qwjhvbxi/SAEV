


function [A,Atimes,ASortInd]=cleanData(AIN,AtimesIN)

if sum(AtimesIN(:,1)~=AtimesIN(:,2))>0
    % impossible trips
    Rmv1=AtimesIN(:,1)>=AtimesIN(:,2);
else
    Rmv1=0;
end

% 0s
Rmv=logical(Rmv1+(AtimesIN(:,1)==0));

% remove wrong trips
AIN(Rmv,:)=[];
AtimesIN(Rmv,:)=[];

% trips sorted by departure time
[Atimes,ASortInd]=sortrows(AtimesIN,1);
A=AIN(ASortInd,:); % trips sorted by departure time


return

% temporary solution for cell files

tripfolder='NYC2018_10wed';
DataFolder=setDataFolder();
A1=A;
Atimes1=Atimes;
for i=1:10
    
    % determine file name
    tripFileLocation=[DataFolder 'trips/' tripfolder '/d' num2str(i) '.mat'];

    A=A1{i};
    Atimes=Atimes1{i};
    save(tripFileLocation,'A','Atimes');
    
end

