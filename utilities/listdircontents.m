%% Files=LISTDIRCONTENTS(folderaddress)
% Print to screen all files in a folder numbered and return the list as
% output. 

function Files=listdircontents(folderaddress)

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