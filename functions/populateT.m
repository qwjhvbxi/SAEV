%% T=populateT(T)
% reconstruct missing coordinates

function T=populateT(T)

addpath('../CarSharingModel/utilities/dijkstra_alg');
N=size(T,1);
for i=1:N
    i
    for j=1:N
        if T(i,j)==0
            T(i,j)=dijkstra(T,i,j);
        end
    end
end

