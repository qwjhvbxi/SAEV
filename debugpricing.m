
if 1

%% testing 

N=[   0.65574 , 0.70605
      0.03571 , 0.53183
      0.84913 , 0.27692
      0.33399 , 0.04617
   ];
% scatter(N(:,1),N(:,2))
% A=[ 0 3 2 0;
%     1 0 9 8;
%     1 1 0 0;
%     3 0 6 0];
% m.v=[5;18;2;9];
A=[ 0 3 2 0;
    3 0 9 8;
    2 9 0 0;
    0 8 0 0];
m.v=[0;0;0;0];
m.a=A;

m.c=round(20*((N(:,1)-N(:,1)').^2+(N(:,2)-N(:,2)').^2));
m.gamma_r=0.1;
m.gamma_p=0.25;
% m.fixedprice=0.5;

[prices,iters,m]=NLPricing(m)

end

% bestp=0.3;


%% checks

c=m.c;
% pr=normalizedprices;

% demand acceptance probability as a function of price
f=@(s) exp(-s.*c)./(exp(-s.*c)+exp(-m.gamma_p*c));

% demand acceptance probability as a function of price and distance
f2=@(s,d) exp(-s.*d)./(exp(-s.*d)+exp(-m.gamma_p*d));

bestprice=(m.gamma_r+m.gamma_p*2)/(4*m.gamma_p);
for j=1:100
    pp=j/100;           % normalized price
    rp=pp*2*m.gamma_p;  % real price
    d=(m.a).*f(rp);       % demand
    revn(j)=sum(sum(d.*c.*rp));
    cost(j)=sum(sum(d.*c*m.gamma_r));
end

d=(m.a).*f(prices);       % demand
revnOpti=sum(sum(d.*c.*prices));
costOpti=sum(sum(d.*c*m.gamma_r));

figure
hold on
plot(revn)
plot(cost)
plot(revn-cost)
line([0,100],[revnOpti-costOpti,revnOpti-costOpti])

%   arrivals      departures     vehicles
% b=sum(full(m.a))-sum(full(m.a'))+m.v';
% m.v=m.v+((b<0).*(-b))';

% function [revn,cost]=calculateresults(m,pp)
% 
% end

return

%%
figure
hold on
p=0:0.01:1;
d=8;
plot(p,f2(p,d).*p)
plot(p,f2(p,d).*p-f2(p,d)*0.1)
plot(p,f2(p,d))


