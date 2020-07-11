%% alternative costs
p2=8;
p3=9;
a=exp(p2)+exp(p3);

%% calculate linearized demand curve
% a= sum of exp of alternatives
% w= distance from center probability for linearization
function [pMin,pMax,b,c]=LinearizedDemandCurve(a,w)
k=0.5; % central probability
s=log(k*a/(1-k)); % point at probability k
b=a*exp(s)./(a+exp(s)).^2; % derivative at s
c=k-b.*s; % intercept of tangent
r=w/b; % left and right linearization distance
pMin=s-r;
pMax=s+r;
end

%% plot
p1=1:0.1:20;
figure
hold on
plot(p1,exp(p1)./(exp(p1)+a))
plot(p1,a*exp(p1)./(a+exp(p1)).^2) % derivative
x=[s-r,s+r];
line(x,b*x+c,'Color','r')





