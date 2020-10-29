%% [X,prices]=RelocationPricing3(m)
% m is a struct with variables: c,v,a,gamma_r,fixedprice
% c is the distance matrix; v is the vehicles at nodes; a is the latent
% demand matrix; gamma_r is the cost of relocation per minute; fixedprice is the
% fixed price for optimizing relocation only (optional).
% 
% X is the relocation matrix, with X(i,j) the number of vehicles moved from
% i to j. The variable is an empty vector when there are no relocations.  
% prices is the price metrix, with prices(i,j) the optimal pricing of a
% trip from i to j. 
%
% see also generalC

function [X,prices,fval]=RelocationPricing3(m)

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
amin=m.amin(:);
amax=m.amax(:);
Da=m.amax-m.amin;
Dp=m.pmax-m.pmin;
Dp(Da==0)=1;


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
Ap=zeros(size(a_ji));
Ap(logical(a_ji))=a.*h;
Ap(logical(a_ij))=Ap(logical(a_ij))-a.*h;

A=[Ar , Ap];
b=m.v+sum(m.a.*s)'-sum(m.a.*s,2);

% constraints on relocation vehicles
Av=[repmat(eye(n),1,n) sparse(n,n^2)];
A=[A;Av];
b=[b;m.v];

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
f=[  c*m.relocationcost ;   -a.*c.*(  s(:) + h*m.relocationcost  )  ];

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

return






