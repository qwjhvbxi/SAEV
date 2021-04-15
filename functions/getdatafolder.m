function DataFolder=getdatafolder()

if exist('functions/DataFolderAddress.mat','file')
    load('functions/DataFolderAddress','DataFolder')
    return
end

if exist('run/datafolderaddress.m','file')
    addpath run
    datafolderaddress
    return
end

DataFolder='data/';
