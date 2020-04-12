%% [Pmat,Resmat]=GENERATEPLOTLINE2(mapscenario,outfieldname,varargin)
% Generate matrices of results from simulations specifying the parameters t
% test. Useful to then plot the values as lines.
% 
% If outfieldname is a char, the Resmat is a vector with the values of the
%   field associated with that name and takes the results from saved files.
%   If results are not found, it returns NaN.
% If outfieldname is a number, the function runs the simulations with the
%   corresponding setting in general5 (second parameter) and returns a
%   struct array of results.  
% if outfieldname is -2, the function deletes the corresponding files
% 
% Warning: variables with text values should be inserted between double
% quotation marks ( " " ) 
% 
% see also GENERAL11, PARAMS

function [Pmat,Resmat]=generateplotline2(mapscenario,outfieldname,varargin)

% default behavior
if nargin<2 || isempty(outfieldname)
    outfieldname=1;
end
if nargin>0  
    if isempty(mapscenario)
        clear mapscenario;
    end
end

numparams=cellfun('length',varargin(2:2:end));
posparams=numparams>1;
whichparm=find(posparams);
varparams=numparams(posparams);
if length(varparams)>2
    Pmat=NaN;
    warning('only a maximum of 2 parameter can vary');
    return
end
if isempty(varparams)
    varparams(1)=1;
    whichparm(1)=0;
end
if length(varparams)==1
    varparams(2)=1;
    whichparm(2)=0;
end

params

N=0;

for k=1:varparams(1)

    for j=1:varparams(2)

        N=N+1;
        
        Pmat{N}=P;

        for i=1:floor(length(varargin)/2)
            
            if whichparm(1)==i
                Pmat{N}.(varargin{(i-1)*2+1})=varargin{(i-1)*2+2}(k);
            elseif whichparm(2)==i
                Pmat{N}.(varargin{(i-1)*2+1})=varargin{(i-1)*2+2}(j);
            else
                Pmat{N}.(varargin{(i-1)*2+1})=varargin{(i-1)*2+2};
            end
        end
    end
end

if ischar(outfieldname)
    for j=1:N
        try
            load(['out/simulations/' DataHash(Pmat{j})],'Res');
            Resmat(j)=Res.(outfieldname);
        catch
            Resmat(j)=NaN;
        end
    end
    Resmat=reshape(Resmat,fliplr(varparams))';
    return
end
    
if outfieldname==-2
	for j=1:N
    	delete(['out/simulations/' DataHash(Pmat{j})],'Res');
	end
    return
end

for j=1:N
    Resmat(j)=general11(Pmat{j},outfieldname,-j); 
end

Resmat=reshape(Resmat,fliplr(varparams))';

end









