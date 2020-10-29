%% [prices,k,m,reloc,revenues]=NLPricing(m)
% Non-linear pricing by OD, with fixed linearization intervals
% 
% m is a struct with variables: c,v,a,gamma_r,gamma_alt,fixedprice
% c is the distance matrix; v is the vehicles at nodes; a is the latent
% demand matrix; gamma_r is the cost of relocation per minute; fixedprice is the
% fixed price for optimizing relocation only (optional).

function [prices,k,m,reloc,revenues]=NLPricing(m)

%% initializations

c=m.c;          % distance matrix
n=size(c,1);    % nodes
c(1:n+1:end)=1; % remove distance between same nodes
m.w=zeros(n,n); % position of price linearization (between -3 and 3)
if isfield(m,'maxiter')
    maxiter=m.maxiter;      % max number of iterations
else
    maxiter=4;      % max number of iterations
end

if ~isfield(m,'altp')
    m.altp=m.gamma_alt.*m.c;
end


%% useful functions

% price at certain probability: given a trip distance d, find the price at
% which the probability of acceptance is x
q=@(d,x) -(log(  (  x.*exp(-m.altp)  )./(1-x)  ))./d;
g=@(a,s) exp(s)./(exp(s)+a);        % value at s
Points=g(exp(0),-3.5:3.5);          % probability linearization intervals (7 intervals, 8 limits)
PL=length(Points)/2;


%% iterations

fprintf('\n iterations: ')

for k=1:maxiter

    fprintf('*');
    
    % probabilities of trip acceptance at the limits
    m.amin=1-Points(PL+1+m.w);
    m.amax=1-Points(PL+m.w);
    
    % limits to price
    m.pmin=q(c,m.amax);
    m.pmax=q(c,m.amin);
    
    % launch pricing optimization for this iteration
    [reloc,prices]=RelocationPricing3(m);

    % detect if there are prices at the edges
    moves=(round(prices,3)==round(m.pmax,3))-(round(prices,3)==round(m.pmin,3));

    ThisRevenue=CalcRevenue(m,prices,reloc);
    revenues(k)=ThisRevenue;
    
    if sum(moves(:)~=0)>0
        m.w=m.w+moves;
    else
        break;
    end

end

% other functions
% f=@(a,s) log(s.*a./(1-s));          % inverse
% d=@(a,s) a*exp(s)./((a+exp(s)).^2); % derivative at s

return


%% debug

m.v=[0 0 10]';
m.c=[0 2 3
     2 0 4
     3 4 0];
m.a=[0 10 0
     5  0 2
     0  0 0];
m.gamma_r=0.1;
m.gamma_alt=0.25;
[prices,k,m]=NLPricing(m)

%% check inputs

j=3;
figure
hold on
line([m.pmin(j),m.pmax(j)],[m.amax(j),m.amin(j)]);
p=0:0.01:1;
plot(p,g(exp(-m.gamma_alt*m.c(j)),-p*m.c(j)),'k:')


