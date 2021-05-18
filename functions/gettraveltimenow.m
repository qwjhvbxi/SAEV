%% [thisT]=GETTRAVELTIMENOW(T,k)
% Return current travel time for scenario with variable travel times.
% T is a struct with fields `hour` and `traveltime`.
% k is the minute in the day (between 0 and 1440)
% 
% See also: mainsim

function [thisT]=gettraveltimenow(T,k)

if isstruct(T)
    k=rem(k,1440);
    traveltimes=[0,T.hour,24.1]*60;
    nextIndex=find(traveltimes>k,1);
    minutePast=traveltimes(nextIndex-1);
    minuteNext=traveltimes(nextIndex);
    hourNextWeight=(k-minutePast)/(minuteNext-minutePast);
    thisT=T(max(1,nextIndex-2)).traveltime*(1-hourNextWeight)+T(min(length(T),nextIndex-1)).traveltime*hourNextWeight;
else
    thisT=T;
    %thisT(1:length(thisT)+1:end)=0;          % no distance between same node
end
