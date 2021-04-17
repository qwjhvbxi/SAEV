%% [Pmat,Resmat]=GENERATEPLOTLINE3(mapscenario,outfieldname,varargin)
% Generate matrices of results from simulations specifying the parameters
% to test. Useful to then plot the values as lines.
% 
% If outfieldname is a char, the Resmat is a vector with the values of the
%   field associated with that name and takes the results from saved files.
%   If results are not found, it returns NaN.
% If outfieldname is a number, the function runs the simulations with the
%   corresponding setting in main (second parameter) and returns a
%   struct array of results.  
% If outfieldname is -2, the function deletes the corresponding files
% 
% After the first 2 inputs, specify the parameters to test with the format
% (...,'FieldName',value,...). You can specify up to two parameters with
% vector values for sensitivity analysis and an unlimited number of
% parameters with scalar values (same for all simulations). 
% 
% Warning: variables with text values should be inserted between double
% quotation marks ( " " ) 
% 
% see also GENERAL11, CPAR

function [Pmat,Resmat]=generateplotline3(mapscenario,outfieldname,varargin)

% default behavior
if nargin<2 || isempty(outfieldname)
    outfieldname=1;
end
if nargin>0  
    if isempty(mapscenario)
        clear mapscenario;
    end
end

DataFolder=getdatafolder();

% generate input matrix of struct
S=getp(mapscenario);
[Pmat,varparams]=changeStruct(S,varargin);

N=numel(Pmat);

if ischar(outfieldname)
    
    FieldName=char(outfieldname);
    PointPos=find(FieldName=='.');

    for j=1:N
        try
            Hash=DataHash(Pmat{j});
            simname=[DataFolder 'out/' Hash '.mat'];
            load(simname,'Res');
            if isempty(PointPos) 
                Resmat(j)=Res.(FieldName);
            else
                Resmat(j)=Res.(FieldName(1:PointPos-1)).(FieldName(PointPos+1:end));
            end
            % Resmat{j}=Res.(outfieldname);
        catch
            Resmat(j)=NaN;
        end
    end
    Resmat=reshape(Resmat,fliplr(varparams))';
    Pmat=reshape(Pmat,fliplr(varparams))';
    return
end
    
if outfieldname==-2
    for j=1:N
        Hash=DataHash(Pmat{j});
        simname=[DataFolder 'out/' Hash '.mat'];
        delete(simname);
    end
    Resmat=NaN;
    return
end

parfor j=1:N
    Resmat(j)=main(Pmat{j},outfieldname,-j); 
end

Resmat=reshape(Resmat,fliplr(varparams))';
Pmat=reshape(Pmat,fliplr(varparams))';

end









