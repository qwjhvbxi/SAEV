%% [Column1,Column2,Resolution]=readexternalfile(FileName,Day)
% return data in each minute for the specified day and the next.

function [Column1,Column2,Resolution]=readexternalfile(FileName,Day,MinuteSample)

fileID=fopen(FileName,'r');
if fileID==-1
    warning('File not found.');
    return;
end
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

if MinuteSample
    % repeat values to find per minute values
    Column1=repelem(T{1}(IntervalAdjusted),Resolution,1);
    Column2=repelem(T{2}(IntervalAdjusted),Resolution,1);
    Resolution=1;
else
    Column1=T{1}(IntervalAdjusted);
    Column2=T{2}(IntervalAdjusted);
end

% replace nan values with zeros
Column1(isnan(Column1))=0;
Column2(isnan(Column2))=0;
