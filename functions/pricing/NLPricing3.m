%% [prices,k,m]=NLPricing(m)
% Non-linear pricing
% 
% m is a struct with variables: c,v,a,gamma_r,gamma_alt,fixedprice
% c is the distance matrix; v is the vehicles at nodes; a is the latent
% demand matrix; gamma_r is the cost of relocation per minute; fixedprice is the
% fixed price for optimizing relocation only (optional).

function [prices,k,m,reloc]=NLPricing3(m)

%% useful functions

% price at certain probability: given a trip distance d, find the price at
% which the probability of acceptance is x
q=@(d,x) -(log(  (  x.*exp(-m.gamma_alt*d)  )./(1-x)  ))./d;
g=@(a,s) exp(s)./(exp(s)+a);        % value at s
Points=g(exp(0),-4.5:4.5);          % probability linearization intervals (7 intervals, 8 limits)
PL=length(Points)/2;

%% initializations

c=m.c;          % distance matrix
n=size(c,1);    % nodes
c(1:n+1:end)=1; % remove distance between same nodes
m.w=zeros(n,n); % position of price linearization (between -3 and 3)
maxIter=5;      % max number of iterations


%% iterations

fprintf('\n iterations: ')

for k=1:maxIter

    fprintf('*');
    
    % probabilities of trip acceptance at the limits
    m.aminlinear=1-Points(PL+1+m.w);
    m.amaxlinear=1-Points(PL+m.w);
    
    % limits to price
    m.pminlinear=q(c,m.amaxlinear);
    m.pmaxlinear=q(c,m.aminlinear);
    
    % alternative
    m.amin=zeros(n,n);
    m.amax=ones(n,n);
    
    m.pmin=(m.pminlinear-m.pmaxlinear+m.amaxlinear.*m.pmaxlinear-m.aminlinear.*m.pminlinear)./(m.amaxlinear-m.aminlinear);
    m.pmax=(m.amaxlinear.*m.pmaxlinear-m.aminlinear.*m.pminlinear)./(m.amaxlinear-m.aminlinear);
    

    %     if 0
    %         % check inputs
    %         j=3;
    %         figure
    %         hold on
    %         line([m.pmin(j),m.pmax(j)],[m.amax(j),m.amin(j)]);
    %         p=0:0.01:1;
    %         plot(p,g(exp(-m.gamma_alt*m.c(j)),-p*m.c(j)),'k:')
    %     end
    
    % launch pricing optimization for this iteration
%     [reloc,prices,~]=RelocationPricing4(m);
    [reloc,prices]=RelocationPricing3(m);

    % detect if there are prices at the edges
    moves=(round(prices,3)>=round(m.pmaxlinear,3))-(round(prices,3)<=round(m.pminlinear,3));

%     if sum((m.w(:).*moves(:))<0)>0
%         m.w
%         moves
%     end
    
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
     3 4 0]*5;
m.a=[0 10 0
     5  0 2
     0  0 0];
m.gamma_r=0.1;
m.gamma_alt=0.25;
m.relocation=0;
[prices,k,m,reloc]=NLPricing3(m)

%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%

n=size(m.c,1);

Fare=m.gamma_alt.*m.c;
Aeff0=( exp(-Fare)./(exp(-Fare)+exp(-m.gamma_alt*m.c))  ).*m.a
Aeff0.*Fare

p1=prices.*m.c;
Aeff=( exp(-p1)./(exp(-p1)+exp(-m.gamma_alt*m.c))  ).*m.a
Aeff.*p1

sum(sum(Aeff0.*(Fare-m.gamma_r*m.c)))
sum(sum(Aeff.*(p1-m.gamma_r*m.c)))

s=m.v+sum(Aeff)'-sum(Aeff,2);
r=sum(reloc)'-sum(reloc,2);

s+r


