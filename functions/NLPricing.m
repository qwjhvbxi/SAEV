%% [prices,k,m]=NLPricing(m)
% Non-linear pricing
% 
% m is a struct with variables: c,v,a,gamma_r,fixedprice
% c is the distance matrix; v is the vehicles at nodes; a is the latent
% demand matrix; gamma_r is the cost of relocation per minute; fixedprice is the
% fixed price for optimizing relocation only (optional).

function [prices,k,m]=NLPricing(m)

%% useful functions

% price at certain probability: given a trip distance d, find the price at
% which the probability of acceptance is x
q=@(d,x) -(log(  (  x.*exp(-m.gamma_alt*d)  )./(1-x)  ))./d;
g=@(a,s) exp(s)./(exp(s)+a);        % value at s
Points=g(exp(0),-3.5:3.5);          % probability linearization intervals (7 intervals, 8 limits)


%% initializations

c=m.c;          % distance matrix
n=size(c,1);    % nodes
c(1:n+1:end)=1; % remove distance between same nodes
m.w=zeros(n,n); % position of price linearization (between -3 and 3)
maxIter=4;      % max number of iterations


%% iterations

fprintf('\n iterations: ')

for k=1:maxIter

    fprintf('*');
    
    % probabilities of trip acceptance at the limits
    m.amin=1-Points(5+m.w);
    m.amax=1-Points(4+m.w);
    
    % limits to price
    m.pmin=q(c,m.amax);
    m.pmax=q(c,m.amin);

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
    [~,prices]=RelocationPricing3(m);

    % detect if there are prices at the edges
    moves=(round(prices,3)==round(m.pmax,3))-(round(prices,3)==round(m.pmin,3));

    if sum(moves(:)~=0)>0
        m.w=m.w+moves;
    else
        break;
    end

end

% other functions
% f=@(a,s) log(s.*a./(1-s));          % inverse
% d=@(a,s) a*exp(s)./((a+exp(s)).^2); % derivative at s



