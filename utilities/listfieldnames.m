function T=listfieldnames(S,T,fieldname)

if nargin<2
    T=table();
end

if isstruct(S)
    fn=fieldnames(S);
    for i=1:length(fn)
        a=S.(fn{i});
        T=listfieldnames(a,T,fn(i));
    end
else
	T=[T;table(fieldname,numel(S))];
end