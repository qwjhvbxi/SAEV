function setupfolders()

if exist('functions/DataFolderAddress.mat','file')
    load('functions/DataFolderAddress.mat','DataFolder');
    if exist('DataFolder','var') && exist(DataFolder,'dir')
        warning('data folder is already specified!')
        return
    end
end

fprintf('\n\n Where to create external folders?');

try 
    DataFolder=uigetdir;
catch
    DataFolder=0;
end

if DataFolder==0
    fprintf('\n\n Specify an absolute address in the form ''c:\\Users\\Documents\\ExtFolder'' \n or a relative address like ''ExtFolder'' ')
    fprintf('\n\n Just press Enter without input to exit setup.\n\n')
    DataFolder=input(' Type in the address: ','s');
end

if ~isempty(DataFolder)
    
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
    mkdir(DataFolder,'par')
    
else
    
    fprintf('\n\n Setup cancelled. \n\n')
    
end