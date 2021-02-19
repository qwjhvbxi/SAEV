function plotwaitingtimes(Atimes,waiting,binsize)

reqtime=ceil(Atimes(:,1)/binsize);
T=max(reqtime);

arrivi=histc(reqtime,1:T);
attese=accumarray(reqtime,waiting);

plot(1:T,attese./arrivi)

end