%% [y]=movingAverage(x,t,binsize)
% calculate moving average of 'x' according to associated time of event
% 't', with bin size 'binsize'. x and t must be the same length. 

function [y]=movingAverage(x,t,binsize)

binsize=round(binsize);
halfbin=floor(binsize/2);
y=zeros(tsim,1);
for k2=1:tsim
    y(k2)=mean(x(logical((t>=k2-halfbin).*(t<=k2+halfbin))),'omitnan');
end
y(isnan(y))=0; % for bins without data, assume 0

end



