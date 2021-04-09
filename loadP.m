function P=loadP(Pname)

addpath functions utilities
DataFolder=setDataFolder();

if nargin<1
    
    Files=displayContents([DataFolder 'par/']);
    
    Pname=input('P var name or number: ','s');
end

PossibleID=str2num(Pname);
if ~isempty(PossibleID) && PossibleID<=length(Files)
    Pfile=Files(PossibleID).name;
else
    Pfile=[Pname '.mat'];
end

Paddr=[DataFolder 'par/' Pfile];

if exist(Paddr,'file')
    load(Paddr,'P');
else 
    warning('File not found.');
    P=[];
end

end



