%% [p1,p2]=ERRBARMEDIAN(x,y,lw,ci,c,s)
% x,y: data
% lw: line width
% ci: if 0, use median and min/max values; if 1, use confidence intervals
% c: line color
% s: line style
% 
% p1,p2 are the plot handles

function [p1,p2]=errbarmedian(x,y,lw,ci,c,s)

if nargin<6
    s='-';
end
if nargin<5
    c='k';
end
if nargin<4
    ci=0;
end
if nargin<3
    lw=1;
end

if ci==0
    m=median(y,'omitnan');
    m1=m-min(y);
    m2=max(y)-m;
end
if ci==1 %|| size(y,2)==1
    [CI,m]=confidenceInterval(y,[0.025,0.975]);
    m1=m-min(CI);
    m2=max(CI)-m;
end

if size(y,2)>1
    p1=errorbar(x,m,m1,m2,'LineWidth',lw,'Color',c,'LineStyle',s);
    p2=[];
else
    hold on
    p1=line(x([1,end])'*ones(1,2),ones(2,1)*[m-m1,m+m2],'LineWidth',lw,'LineStyle',':','Color',c); 
    p2=line(x([1,end])'',ones(2,1)*m,'LineWidth',lw,'Color',c,'LineStyle',s); 
end

end