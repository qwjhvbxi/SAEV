%% alternative costs
p2=-20;
p3=-12;
a=exp(p2)+exp(p3);
w=0.3;

%% calculate linearized demand curve
% note: non deve essere alla meta'! il punto deve essere determinato dalla tariffa base!!

% a= sum of exp of alternatives
% w= distance from center probability for linearization
% function [pMin,pMax,b,c]=LinearizedDemandCurve(a,w)

q=@(k) log(k*a/(1-k)); % point at probability k
s=q(0.5); % central probability

b=a*exp(s)./(a+exp(s)).^2;  % derivative at s
b=(0.9-0.1)/(q(0.9)-q(0.1));      % derivative between 0.1 and 0.9

c=k-b.*s; % intercept of tangent
r=w/b; % left and right linearization distance
pMin=s-r;
pMax=s+r;
% end

%% plot
p1=-20:0.1:0;
figure
hold on
plot(p1,exp(p1)./(exp(p1)+a))
plot(p1,a*exp(p1)./(a+exp(p1)).^2) % derivative
x=[s-r,s+r];
line(x,b*x+c,'Color','r')
xlabel('-price')
ylabel('demand')





