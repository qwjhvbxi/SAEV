%% [prices,k,m]=NLPricing5(m)
% Non-linear pricing with continuous approximation
% 
% m is a struct with variables: c,v,a,relocationcost,gamma_alt,fixedprice
% c is the distance matrix; v is the vehicles at nodes; a is the latent
% demand matrix; relocationcost is the cost of relocation per minute; fixedprice is the
% fixed price for optimizing relocation only (optional).

function [prices,k,m,reloc,revenues]=NLPricing5(m)


%% initializations

m.c=max(1,m.c);
n=size(m.c,1);    % nodes
if isfield(m,'maxiter')
    maxiter=m.maxiter;      % max number of iterations
else
    maxiter=5;      % max number of iterations
end

if ~isfield(m,'altp')
    m.altp=m.gamma_alt.*m.c;
end

prices=m.altp./m.c;

% revenues=zeros(maxiter+1,1);
% ThisRevenue=CalcRevenue(m,prices);
% revenues(1)=ThisRevenue;

%% useful functions

% price at certain probability: given a trip distance d, find the price at
% which the probability of acceptance is x
q=@(d,x) -(log(  (  x.*exp(-m.altp)  )./(1-x)  ))./d;
g=@(a,s) exp(s)./(exp(s)+a);        % value at s
d=@(a,z,c) -(a.*c.*exp(z.*c))./((1+a.*exp(z.*c)).^2); % derivative at s
Points=g(exp(0),-4.5:4.5);          % probability linearization intervals (7 intervals, 8 limits)
PL=length(Points)/2;


%% iterations

delta=zeros(n^2,maxiter);
priceshist=zeros(n^2,maxiter);

fprintf('\n iterations: ')

for k=1:maxiter

    % coefficients of y=Dx+C
    D=d(exp(-m.altp),prices,m.c);
    C=g(exp(-m.altp),-prices.*m.c)-D.*prices;
    
    % alternative
    
    p1=(1-C)./D;
    p0=-C./D;
    
    % m.pmin=p1;
    % m.pmin=prices-(prices-p1)/k;
    m.pmin=prices-0.5./m.c; 

    if isfield(m,'fixedfleet') && m.fixedfleet
        m.pmax=p0;
    else
        m.pmax=prices+0.5./m.c;
    end
    
    m.amin=D.*m.pmax+C;%g(exp(-m.altp),-m.pmax.*m.c);
    m.amax=D.*m.pmin+C;%g(exp(-m.altp),-m.pmin.*m.c);
    
    Empty=(m.a==0);
    m.amin(Empty)=0;
    m.amax(Empty)=0;
    m.pmin(Empty)=0;
    m.pmax(Empty)=0;

    % launch pricing optimization for this iteration
    [reloc,newprices,fval]=RelocationPricing9(m);
    
    delta(:,k)=(newprices(:)-prices(:));
    priceshist(:,k)=newprices(:);
    
    prices=newprices;
    
    fprintf('*');
    % fprintf('\n %f',mean(delta(:,k)));
    % fprintf('\n %f',fval);
    % ThisRevenue=CalcRevenue(m,prices,reloc);
    % fprintf('\n %f',ThisRevenue);
    % revenues(k)=ThisRevenue;
    
end

% delta1=delta(sum(delta,2)>0,:);
% figure
% plot(mean(delta))
% figure
% plot(delta1')

end

function ThisRevenue=CalcRevenue(m,prices,reloc)
    p1=prices.*m.c;
    Aeff=( exp(-p1)./(exp(-p1)+exp(-m.altp))  ).*m.a;
    if ~isempty(reloc)
        relocation=sum(sum(reloc.*m.c))*m.relocationcost;
    else 
        relocation=0;
    end
    ThisRevenue=full(sum(sum(Aeff.*(p1-m.relocationcost*m.c))))-relocation;
end


function debug1

m.v=[0 0 0]';
m.c=[0 2 3
     2 0 4
     3 4 0]*5;
m.a=[0 10 0
     5  0 0
     0  2 0];
m.relocationcost=0.1;
m.gamma_alt=0.25;
m.relocation=0;
[prices,k,m,reloc]=NLPricing4(m)

%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%

n=size(m.c,1);

Fare=m.gamma_alt.*m.c;
Aeff0=( exp(-Fare)./(exp(-Fare)+exp(-m.gamma_alt*m.c))  ).*m.a
Aeff0.*Fare

p1=prices.*m.c;
Aeff=( exp(-p1)./(exp(-p1)+exp(-m.gamma_alt*m.c))  ).*m.a
Aeff.*p1

sum(sum(Aeff0.*(Fare-m.relocationcost*m.c)))
sum(sum(Aeff.*(p1-m.relocationcost*m.c)))

s=m.v+sum(Aeff)'-sum(Aeff,2);
r=sum(reloc)'-sum(reloc,2);

s+r


% check inputs
j=3;
figure
hold on
line([m.pmin(j),m.pmax(j)],[m.amax(j),m.amin(j)]);
p=0:0.01:1;
plot(p,g(exp(-m.altp(j)),-p*m.c(j)),'k:')

end

