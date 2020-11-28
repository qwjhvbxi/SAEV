%% [elep,co2]=ReadExtFile(FileName,Day)
% return electricity price and carbon intensity in each minute for the
% specified day and the next.

function [elep,co2]=ReadExtFile(FileName,Day)

fileID=fopen(FileName,'r');
H=textscan(fileID,'%q %f',1,'Delimiter',','); % retrieve resolution info
T=textscan(fileID,'%f %f','Delimiter',',','HeaderLines',2); % retrieve values
fclose(fileID);

% resolution in minutes
Resolution=H{2};

% how many data points in each day
PointsPerDay=1440/Resolution;

% total days covered
TotalDays=floor(length(T{1})/PointsPerDay);

% interval of interest
Interval=PointsPerDay*(Day-1)+1:PointsPerDay*(Day+1);

% adjust for end of the dataset (the next day from last day is back to the first)
IntervalAdjusted=rem(Interval-1,TotalDays*PointsPerDay)+1;

% repeat values to find per minute values
elep=repelem(T{1}(IntervalAdjusted),Resolution,1);
co2=repelem(T{2}(IntervalAdjusted),Resolution,1);

% replace nan values with zeros
co2(isnan(co2))=0;


