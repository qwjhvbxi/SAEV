
function clusters=assignclusters(C,centroids)

if size(centroids,2)==1
    D=generatedistancemat(C);
    [~,clusters]=min(D(:,centroids),[],2);
elseif size(centroids,2)==2
    D=generatedistancemat(C,centroids);
    [~,clusters]=min(D,[],2);
end

end