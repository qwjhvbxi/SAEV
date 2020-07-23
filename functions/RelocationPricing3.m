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

% calculations 
n=size(m.c,1);    % nodes
a=reshape(m.a,n^2,1);
c=reshape(m.c,n^2,1);
pmin=reshape(m.pmin,n^2,1);
pmax=reshape(m.pmax,n^2,1);
amin=reshape(m.amin,n^2,1);
amax=reshape(m.amax,n^2,1);

% demand is in the form d=a*(s-h*p)
s=m.amax+(m.amax-m.amin)./(m.pmax-m.pmin).*m.pmin; % static term of demand
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
lb=[zeros(n^2,1);   pmin];
ub=[repmat(m.v,n,1);pmax];
ub(1:n+1:n^2)=0; % no relocation in same node

if isfield(m,'fixedprice') && ~isempty(m.fixedprice)

    % bounds for fixed pricing
    lb(n^2+1:2*n^2)=m.fixedprice;
    ub(n^2+1:2*n^2)=m.fixedprice;
    
end

% cost function
H=2*diag([zeros(n^2,1);a.*c.*h*m.gamma_p]);
f=[c*m.gamma_r;-a.*c.*s_v*m.gamma_p-a.*c.*h*m.gamma_r];
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






%% testing

% N=[   0.65574 , 0.70605
%       0.03571 , 0.03183
%       0.84913 , 0.27692
%       0.93399 , 0.04617
%       0.67874 , 0.09713
%       0.75774 , 0.82346
%       0.74313 , 0.69483
%       0.39223 , 0.31710
%       0.65548 , 0.95022
%       0.17119 , 0.03444
%       ];
% n=length(N);    % nodes
% m.a_ts=rand(n,n)*10;
% m.a_to=m.a_ts+rand(n,n)*10;
% m.v=rand(n,1)*40;

N=round([   0.65574 , 0.70605
      0.03571 , 0.03183
      0.84913 , 0.27692
      0.93399 , 0.04617
      ]*30);
A=[ 0 3 2 0;
    1 0 9 8;
    1 1 0 0;
    3 0 6 0];
m.a=A;
m.v=[5;18;2;9];

m.c=(N(:,1)-N(:,1)').^2+(N(:,2)-N(:,2)').^2;
m.gamma_r=0.1;
m.gamma_p=0.5;
% m.fixedprice=0.5;

[relocations,prices]=RelocationPricing3(m)


