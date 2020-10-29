function setupFolders()

fprintf('\n\n Where to create external folders?');

try 
    DataFolder=uigetdir;
catch
    DataFolder=0;
end

if DataFolder==0
    fprintf('Specify an absolute address in the form ''c:\\Users\\Documents\\'' \n or a relative address \n\n')
    DataFolder=input(' Type in the address: ','s');
    
end

if ~strcmp(DataFolder(end),'/') && ~strcmp(DataFolder(end),'\')
    DataFolder=[DataFolder '\'];
end

save('functions/DataFolderAddress','DataFolder')

if ~exist(DataFolder,'dir')
    mkdir(DataFolder);
end

mkdir(DataFolder,'trips')
mkdir(DataFolder,'scenarios')
mkdir(DataFolder,'out_saev')
mkdir(DataFolder,'figures')
mkdir(DataFolder,'eleprices')
mkdir(DataFolder,'grid')