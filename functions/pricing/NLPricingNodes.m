%% [prices,k,m]=NLPricing(m)
% Non-linear pricing at nodes
% 
% m is a struct with variables: c,v,a,relocationcost,gamma_alt,fixedprice
% c is the distance matrix; v is the vehicles at nodes; a is the latent
% demand matrix; relocationcost is the cost of relocation per minute; fixedprice is the
% fixed price for optimizing relocation only (optional).

function [prices,k,m,reloc]=NLPricingNodes(m)

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

% price at certain probability: given a trip distance d, find the
% additional surcharge at which the probability of acceptance is x
q=@(d,x) -log(  (  x.*exp(-m.altp)  )./(1-x)  ) - m.gamma_d(d).*d;
s=@(d,S) exp(-m.gamma_d(d).*d-S)./(exp(-m.gamma_d(d).*d-S)+exp(-m.altp));
g=@(a,s) exp(s)./(exp(s)+a);        % value at s
u=@(a,x) log((x.*a)./(1-x));        % value at s


%% iterations

% probability at 0 surcharge
Surcharges=zeros(n);
m.pvec=zeros(1,n*2);

fprintf('\n iterations: ')

for k=1:maxiter

    fprintf('*');
    
    % limits to price
    m.pmin=Surcharges-0.5;
    m.pmax=Surcharges+0.5;
    
    % probabilities of trip acceptance at the limits
    m.amin=s(c,m.pmax);
    m.amax=s(c,m.pmin);
    m.amin(1:n+1:end)=0.5;
    m.amax(1:n+1:end)=0.5;

    [reloc,prices]=RelocationPricing6(m);
    
    if ~isempty(prices)
        Surcharges=prices(1:n)+prices(n+1:n*2)';
        m.pvec=prices';
    else
        prices=m.pvec';
        break
    end
    
end

return






%% debug

m.v=[10 0 0]';
m.c=[0 15 29
     15 0 14
     29 14 0];
m.a=[0 10 0
     5  0 0
     0  2 0];
% m.a=[0 10 0
%      10 0 0
%      0  0 0];
m.relocationcost=0.1;
m.gamma_alt=0.25;
m.gamma_d=bestp((1:50)',m.relocationcost,m.gamma_alt);
[prices,k,m,reloc]=NLPricing2(m)

%%

n=size(m.c,1);
Surcharges=prices(1:n)+prices(n+1:2*n)';

Fare=(m.gamma_d(m.c+eye(n))).*m.c;
Aeff0=( exp(-Fare)./(exp(-Fare)+exp(-m.gamma_alt*m.c))  ).*m.a
Aeff0.*Fare

p1=Fare+Surcharges;
Aeff=( exp(-p1)./(exp(-p1)+exp(-m.gamma_alt*m.c))  ).*m.a
Aeff.*p1

sum(sum(Aeff0.*Fare-m.relocationcost*m.c))
sum(sum(Aeff.*p1-m.relocationcost*m.c))

s=m.v+sum(Aeff)'-sum(Aeff,2);
r=sum(reloc)'-sum(reloc,2);

s+r

%% check inputs

j=3;
figure
hold on
line([m.pmin(j),m.pmax(j)],[m.amax(j),m.amin(j)]);
p=-2:0.1:2;
plot(p,  g(  exp(-m.gamma_alt*m.c(j))  ,   -m.gamma_d(m.c(j))*m.c(j)-p)  ,'k:')

%% other

[MinSorted,Ind1]=sort(m.pmin(:));
plot([MinSorted,m.pmax(Ind1)])
