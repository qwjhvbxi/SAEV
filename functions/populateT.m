%% T=populateT(T)
% reconstruct missing coordinates (==0) with Dijkstra algorithm

function T=populateT(T)

addpath('utilities/dijkstra_alg');
N=size(T,1);
for i=1:N
    i
    for j=1:N
        if T(i,j)==0
            T(i,j)=dijkstra(T,i,j);
        end
    end
end

