


function [prices,k,m]=NLPricing(m)

% utilities
c=m.c;
n=size(c,1);  % nodes
c(1:n+1:end)=1;
g=@(a,s) exp(s)./(exp(s)+a);        % value at s
f=@(a,s) log(s.*a./(1-s));          % inverse
d=@(a,s) a*exp(s)./((a+exp(s)).^2); % derivative at s
Points=g(exp(0),-3.5:3.5);          % probability linearization intervals (7 intervals, 8 limits)

% price at certain probability: given a trip distance d, find the price at
% which the probability of acceptance is x
q=@(d,x) -(log(x.*exp(-m.gamma_p*d)./(1-x)))./d;

% initializations
% a=exp(bp*c); % exp of benefit of alternative modes for each node pair
m.w=zeros(n,n);    % position of linearization (between -3 and 3)
u=0.5;
maxIter=2;

for k=1:maxIter

    k
    
    % probabilities of trip acceptance at the limits
    m.amin=1-Points(5+m.w);
    m.amax=1-Points(4+m.w);
    
    % limits to price
    m.pmin=q(c,m.amax);
    m.pmax=q(c,m.amin);

    [reloc,prices]=RelocationPricing3(m);

    movedown=(round(prices,2)==round(m.pmin,2));
    moveup=(round(prices,2)==round(m.pmax,2));

    if sum(movedown(:))+sum(moveup(:))>0
        m.w=m.w-movedown;
        m.w=m.w+moveup;
    else
        break;
    end

end

% reloc

% fmin=f(a,2*p2*pmin0.*c);
% fmax=f(a,2*p2*pmax0.*c);
% 
% fmin=f(a,p2*c+w-u);
% fmax=f(a,p2*c+w+u);
% pmax=1-fmin;
% pmin=1-fmax;
% pmax(1:n+1:end)=0.5;
% pmin(1:n+1:end)=0.5;
% 
% m.pmin=
% m.pmax=


return



