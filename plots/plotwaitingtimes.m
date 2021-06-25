function plotwaitingtimes(Atimes,waiting,binsize)

if length(Atimes)<length(waiting)
    % cumulativeTripArrivals
    
    T=ceil(1440/binsize);
    
    Atimes=Atimes(1:binsize:end);
    arrivi=Atimes(2:T+1)-Atimes(1:T);
    
    attese=zeros(T,1);
    for i=1:T
        
        attese(i)=sum(waiting(Atimes(i)+1:min(length(waiting),Atimes(i+1))));
        
    end
    
else

    % request time bin
    reqtime=ceil(Atimes(:,1)/binsize);

    % total sim time bins
    T=max(reqtime);

    % number of arrivals for each time bin
    arrivi=histc(reqtime,1:T);

    % waiting time in each bin
    attese=accumarray(reqtime,waiting);

end

plot(1:T,attese./arrivi)

end