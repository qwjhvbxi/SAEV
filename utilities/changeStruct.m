%% [Pmat,varparams]=changeStruct(S,Inputs)
% Return a matrix of struct with modified parameters. The parameters to
% change are specified in the format (...,'FieldName',value,...). You can
% specify up to two parameters with vector values and an unlimited number
% of parameters with scalar values (same for all simulations).  
% 
% Variables with text values should be inserted between double
% quotation marks ( " " ).
%   example: (...,'FieldName',"text",...)

function [Pmat,varparams]=changeStruct(S,varargin)

if numel(varargin)==1
    Inputs=varargin{1};
else
    Inputs=varargin;
end

% numparams=cellfun('length',Inputs(2:2:end)); % length of each set of parameter values
numparams=zeros(1,floor(length(Inputs)/2));
for i=1:length(numparams)
    numparams(i)=length(Inputs{2*i});
end
posparams=numparams>1;      % set with more than one values (sensitivity)
whichparm=find(posparams); 
varparams=numparams(posparams); % number of sets with more than 1 value
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

N=0;

for k=1:varparams(1)

    for j=1:varparams(2)

        N=N+1;
        
        Pmat{N}=S;

        for i=1:floor(length(Inputs)/2)
            
            if whichparm(1)==i
                ThisValue=Inputs{(i-1)*2+2}(k);
            elseif whichparm(2)==i
                ThisValue=Inputs{(i-1)*2+2}(j);
            else
                ThisValue=Inputs{(i-1)*2+2};
            end
            
            FieldNameTot=char(Inputs{(i-1)*2+1});
            FieldNames=split(FieldNameTot,',');
            
            if isstring(ThisValue)
                ThisValue=char(ThisValue);
            end
            
            for VarL=1:length(FieldNames)
            
                FieldName=FieldNames{VarL};
                PointPos=find(FieldName=='.');
            
                % todo: check if ThisValue is string
                if isempty(PointPos) 
                    Pmat{N}.(FieldName)=ThisValue;
                else
                    Pmat{N}.(FieldName(1:PointPos-1)).(FieldName(PointPos+1:end))=ThisValue;
                end
            end
            
        end
    end
end