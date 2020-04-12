function [CI,M]=confidenceInterval(x,confint)

if nargin<2
    confint=[0.05 0.95]; % 90% confidence interval
end

SEM=std(x,'omitnan')./sqrt(size(x,1));% Standard Error
ts=tinv(confint,size(x,1)-1);         % T-Score
CI=mean(x,'omitnan') + ts'*SEM;       % Confidence Intervals
M=mean(x,'omitnan');

end