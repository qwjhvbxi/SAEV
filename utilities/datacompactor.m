function b=datacompactor(a)

if sum(a(:)==0)>0.9*length(a(:))
    b=sparse(a);
else
    b=single(a);
end

end