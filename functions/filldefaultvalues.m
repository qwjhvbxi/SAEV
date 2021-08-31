function Q=filldefaultvalues(D,S)

Q=D;

if isstruct(S)
    Sfields=fieldnames(S);
    for i=1:length(Sfields)
        Q.(Sfields{i})=S.(Sfields{i});
    end
end

end