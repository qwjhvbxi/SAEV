


function [A,Atimes,ASortInd]=cleantripdata(AIN,AtimesIN,RmSameNode)

Rmv1=0;

if size(AtimesIN,2)>1
    
    if sum(AtimesIN(:,1)~=AtimesIN(:,2))>0
        % impossible trips
        Rmv1=AtimesIN(:,1)>=AtimesIN(:,2);
    end

end

% 0s
Rmv=logical(Rmv1+(AtimesIN(:,1)==0));

% remove wrong trips
AIN(Rmv,:)=[];
AtimesIN(Rmv,:)=[];


if nargin>2 && RmSameNode==true

    % remove trips inside same node
    SameNode=logical(AIN(:,1)==AIN(:,2));
    AIN(SameNode,:)=[];
    AtimesIN(SameNode,:)=[];

end

% trips sorted by departure time
[Atimes,ASortInd]=sortrows(AtimesIN,1);
A=AIN(ASortInd,:); % trips sorted by departure time


return

% temporary solution for cell files

tripfolder='NYC2018_10wed';
DataFolder=getdatafolder()
A1=A;
Atimes1=Atimes;
for i=1:10
    
    % determine file name
    tripFileLocation=[DataFolder 'trips/' tripfolder '/d' num2str(i) '.mat'];

    A=A1{i};
    Atimes=Atimes1{i};
    save(tripFileLocation,'A','Atimes');
    
end

