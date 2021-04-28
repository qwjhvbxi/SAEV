function b=datacompactor(a)

if nnz(a)<0.1*numel(a)
    b=sparse(a);
else
    b=single(a);
end

end