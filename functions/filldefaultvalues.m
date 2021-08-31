function Q=filldefaultvalues(D,S)

if isstruct(D) && isstruct(S)
    Q=D;
    Sfields=fieldnames(S);
    for i=1:length(Sfields)
        if isfield(Q,Sfields{i})
            Q.(Sfields{i})=filldefaultvalues(Q.(Sfields{i}),S.(Sfields{i}));
        else
            Q.(Sfields{i})=S.(Sfields{i});
        end
    end
else
    Q=S;
end

end