%% [X,prices]=RelocationPricing3(m)
% m is a struct with variables: c,v,a,gamma_r,gamma_p,fixedprice
% c is the distance matrix; v is the vehicles at nodes; a is the latent
% demand matrix; gamma_r is the cost of relocation
% per minute; gamma_p is the base tariff per passenger; fixedprice is the
% fixed price for optimizing relocation only (optional).
% 
% X is the relocation matrix, with X(i,j) the number of vehicles moved from
% i to j. The variable is an empty vector when there are no relocations.  
% prices is the price metrix, with prices(i,j) the optimal pricing of a
% trip from i to j. 
%
% see also generalC

function [X,prices]=RelocationPricing3(m)

% initializations
n=size(m.c,1);    % nodes
a=reshape(m.a,n^2,1);
c=reshape(m.c,n^2,1);
pmin=reshape(m.pmin,n^2,1);
pmax=reshape(m.pmax,n^2,1);
amin=reshape(m.amin,n^2,1);
amax=reshape(m.amax,n^2,1);

% demand is in the form d=a*(s-h*p)
s=m.amax+((m.amax-m.amin)./(m.pmax-m.pmin)).*m.pmin; % static term of demand
s_v=reshape(s,n^2,1);
h=(amax-amin)./(pmax-pmin); % p multiplier for demand

% constraint on relocation
Ar=[repmat(eye(n),1,n)-kron(eye(n),ones(1,n))];

a_ji=kron(eye(n),ones(1,n));
a_ij=repmat(eye(n),1,n);
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
lb=[zeros(n^2,1);    min(1,max(0,pmin))];
ub=[repmat(m.v,n,1); max(0,min(1,pmax))];
ub(1:n+1:n^2)=0; % no relocation in same node

if isfield(m,'fixedprice') && ~isempty(m.fixedprice)

    % bounds for fixed pricing
    lb(n^2+1:2*n^2)=m.fixedprice;
    ub(n^2+1:2*n^2)=m.fixedprice;
    
end

% cost function
H=2*diag([zeros(n^2,1);a.*c.*h*m.gamma_p*2]);
f=[  c*m.gamma_r ; ...
    -a.*c.*(s_v*m.gamma_p*2+h*m.gamma_r)  ];
% f=[c_v*m.gamma_r;-a_to_v*m.gamma_p];

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






