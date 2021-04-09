%% T=populateT(T)
% reconstruct missing coordinates (==0) with Dijkstra algorithm

function newT=populateT(T)

addpath('utilities/dijkstra_alg');
N=size(T,1);
newT=T;
parfor i=1:N
    i
    for j=1:N
        if i~=j
            if T(i,j)==0
                newT(i,j)=dijkstra(T,i,j);
            end
        end
    end
end

