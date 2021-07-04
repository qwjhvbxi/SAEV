%% [prices,k,m,reloc]=NLPRICINGOD(m)
% Non-linear pricing with continuous approximation
% 
% m is a struct with variables:
% c,v,a,altp,movingcostkm,fixedprice,maxiter
% c is the distance matrix; v is the vehicles at nodes; a is the latent
% demand matrix; movingcostkm is the cost of relocation per km; fixedprice is the
% fixed price for optimizing relocation only (optional).

function [prices,k,m,reloc,revenues]=nlpricingod(m)


%% initializations

m.c=max(0.1,m.c);   % distance matrix (km)
n=size(m.c,1);      % number of nodes
prices=m.altp./m.c; % initial prices

% revenues=zeros(maxiter+1,1);
% ThisRevenue=CalcRevenue(m,prices);
% revenues(1)=ThisRevenue;


%% useful functions

% value at s
g=@(a,s) exp(s)./(exp(s)+a);       

% derivative at s
d=@(a,z,c) -(a.*c.*exp(z.*c))./((1+a.*exp(z.*c)).^2); 


%% iterations

delta=zeros(n^2,m.maxiter);
priceshist=zeros(n^2,m.maxiter);

fprintf('\n iterations: ')

for k=1:m.maxiter

    % coefficients of y=Dx+C
    D=d(exp(-m.altp),prices,m.c);
    C=g(exp(-m.altp),-prices.*m.c)-D.*prices;
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


function [X,prices,fval]=RelocationPricing9(m)
%% [X,prices]=RelocationPricing9(m)
% m is a struct with variables: c,v,a,gamma_r,fixedprice
% c is the distance matrix; v is the vehicles at nodes; a is the latent
% demand matrix; gamma_r is the cost of relocation per minute; fixedprice is the
% fixed price for optimizing relocation only (optional).
% 
% X is the relocation matrix, with X(i,j) the number of vehicles moved from
% i to j. The variable is an empty vector when there are no relocations.  
% prices is the price metrix, with prices(i,j) the optimal pricing of a
% trip from i to j. 

%% input check

InputCheck=((m.pmax>=m.pmin).*(m.amax>=m.amin));
if prod(InputCheck(:))==0
    warning('impossible')
    return
end


%% initializations

n=size(m.c,1);    % nodes
a=m.a(:);
c=m.c(:);
pmin=m.pmin(:);
pmax=m.pmax(:);
Da=m.amax-m.amin;
Dp=m.pmax-m.pmin;
Dp(Da==0)=1;
V=sum(m.v); % fleet size


%% demand function
% demand is in the form d=a*(s-h*p)

s=m.amax+(Da./Dp).*m.pmin; % static term of demand
h0=Da./Dp; % p multiplier for demand
h=h0(:);

% checks
% m.a.*(s-h.*m.pmin)<=m.a;
% m.a.*(s-h.*m.pmax)>=0;
% j=3;
% p=pmin(j):0.01:pmax(j);
% dapprox=a(j)*(s(j)-h(j)*p)
% dreal=a(j)*exp(-p*m.c(j))./(exp(-p*m.c(j))+exp(-0.25*m.c(j)))

if isfield(m,'relocation')
    Relocation=m.relocation;
else
    Relocation=1;
end

%% constraints

% constraint on relocation
Ar=[repmat(eye(n),1,n)-kron(eye(n),ones(1,n))];

a_ji=kron(eye(n),ones(1,n)); % sum_j (a_ji)
a_ij=repmat(eye(n),1,n);     % sum_j (a_ij)
Ap=zeros(size(a_ij));
Ap(logical(a_ji))=a.*h;
Ap(logical(a_ij))=Ap(logical(a_ij))-a.*h;

A=[Ar , Ap];
b=sum(m.a.*s)'-sum(m.a.*s,2);

if isfield(m,'fixedfleet') && m.fixedfleet

    % constraints on fleet size
    Av=[ones(1,n^2) , -(a.*h)'];
    bv=-sum(a.*s(:))+V;
    A=[A;Av];
    b=[b;bv];
    
end

% bounds
lb=[zeros(n^2,1);    pmin];
ub=[repmat(m.v,n,1)*Relocation; pmax];
ub(1:n+1:n^2)=0; % no relocation in same node

if isfield(m,'fixedprice') && ~isempty(m.fixedprice)

    % bounds for fixed pricing
    lb(n^2+1:2*n^2)=m.fixedprice;
    ub(n^2+1:2*n^2)=m.fixedprice;
    
end

% cost function
H=2*diag([zeros(n^2,1);a.*c.*h]);
f=[  c*m.movingcostkm ;   -a.*c.*(  s(:) + h*m.movingcostkm  )  ];

options=optimoptions('quadprog','display','none');

% optimization
[X0,fval]=quadprog(H,f,A,b,[],[],lb,ub,[],options);

% results (shift from convention [d,o] to [o,d])
X1=round(X0,3);
X=reshape(X1(1:n^2),n,n);
prices=reshape(X1(n^2+1:end),n,n);

% there are no relocation actions
if sum(X(:))==0
    X=[];
end

end


function ThisRevenue=CalcRevenue(m,prices,reloc)
    p1=prices.*m.c;
    Aeff=( exp(-p1)./(exp(-p1)+exp(-m.altp))  ).*m.a;
    if ~isempty(reloc)
        relocation=sum(sum(reloc.*m.c))*m.movingcostkm;
    else 
        relocation=0;
    end
    ThisRevenue=full(sum(sum(Aeff.*(p1-m.movingcostkm*m.c))))-relocation;
end


function debug1

m.v=[0 0 0]';
m.c=[0 2 3
     2 0 4
     3 4 0]*5;
m.a=[0 10 0
     5  0 0
     0  2 0];
m.movingcostkm=0.1;
m.altp=0.25*m.c;
m.relocation=0;
[prices,k,m,reloc]=nlpricingod(m)

%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%===%

n=size(m.c,1);

Aeff0=( exp(-m.altp)./(exp(-m.altp)+exp(-m.altp))  ).*m.a
Aeff0.*m.altp

p1=prices.*m.c;
Aeff=( exp(-p1)./(exp(-p1)+exp(-m.altp))  ).*m.a
Aeff.*p1

sum(sum(Aeff0.*(m.altp-m.movingcostkm*m.c)))
sum(sum(Aeff.*(p1-m.movingcostkm*m.c)))

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

