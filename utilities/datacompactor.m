%% b=DATACOMPACTOR(a)
% Turns input vector 'a' into sparse array or single precision to save
% space.

function b=datacompactor(a)

if nnz(a)<0.1*numel(a)
    b=sparse(a);
else
    b=single(a);
end

end