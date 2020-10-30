function P=loadP(Pname)

addpath functions
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




function Files=displayContents(folderaddress)

DirFiles=dir(folderaddress);
q=0;

fprintf('Files in folder: \n\n')

for i=1:length(DirFiles)
    if ~strcmp(DirFiles(i).name,'.') && ~strcmp(DirFiles(i).name,'..') && DirFiles(i).isdir==0
        q=q+1;
        Files(q)=DirFiles(i);
        fprintf('%d. %s \n',q,DirFiles(i).name(1:end-4));
    end
end

if q==0
    fprintf('Empty folder. \n\n')
    Files=[];
else
    fprintf('\n')
end
end