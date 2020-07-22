%% [X,prices]=RelocationPricing(m)
% m is a struct with variables: c,v,a_ts,a_to,gamma_r,gamma_p,fixedprice
% c is the distance matrix; v is the vehicles at nodes; a_ts is the latent
% demand matrix during first ts minutes; a_to is the latent demand matrix
% during first to minutes, with to>ts; gamma_r is the cost of relocation
% per minute; gamma_p is the base tariff per passenger; fixedprice is the
% fixed price for optimizing relocation only (optional).
% 
% X is the relocation matrix, with X(i,j) the number of vehicles moved from
% i to j. The variable is an empty vector when there are no relocations.  
% prices is the price metrix, with prices(i,j) the optimal pricing of a
% trip from i to j. 
%
% see also generalC

function [X,prices]=RelocationPricing(m)

% calculations (shift from convention [o,d] to [d,o])
n=size(m.c,1);    % nodes
a_ts_v=reshape(m.a_ts',n^2,1);
a_to_v=reshape(m.a_to',n^2,1);
c_v=reshape(m.c',n^2,1);

% constraint on relocation
a0to=kron(eye(n),ones(1,n));
a0ts=repmat(eye(n),1,n);
a2=zeros(size(a0to));
a2(logical(a0to))=-a_to_v;
a2(logical(a0ts))=a2(logical(a0ts))+a_ts_v;
A=[kron(eye(n),ones(1,n))-repmat(eye(n),1,n)  ,  a2];
b=m.v+sum(m.a_ts,2)-sum(m.a_to,1)';

% bounds
lb=zeros(n^2*2,1);
ub=[repelem(m.v,n,1);ones(n^2,1)];
ub(1:n+1:n^2)=0; % no relocation in same node

if isfield(m,'fixedprice') && ~isempty(m.fixedprice)

    % bounds for fixed pricing
    lb(n^2+1:2*n^2)=m.fixedprice;
    ub(n^2+1:2*n^2)=m.fixedprice;
    
end

% cost function
H=2*diag([zeros(n^2,1);a_to_v.*c_v*m.gamma_p]);
f=[c_v*m.gamma_r;-a_to_v.*c_v*m.gamma_p-a_to_v.*c_v*m.gamma_r];
% f=[c_v*m.gamma_r;-a_to_v*m.gamma_p];

options=optimoptions('quadprog','display','none');

% optimization
X0=quadprog(H,f,A,b,[],[],lb,ub,[],options);

% results (shift from convention [d,o] to [o,d])
X1=round(X0,3);
X=reshape(X1(1:n^2),n,n)';
prices=reshape(X1(n^2+1:end),n,n)';

% there are no relocation actions
if sum(X(:))==0
    X=[];
end

return






%% testing
 
N=[   0.65574 , 0.70605
      0.03571 , 0.03183
      0.84913 , 0.27692
      0.93399 , 0.04617
      0.67874 , 0.09713
      0.75774 , 0.82346
      0.74313 , 0.69483
      0.39223 , 0.31710
      0.65548 , 0.95022
      0.17119 , 0.03444
      ];
n=length(N);    % nodes

m.c=(N(:,1)-N(:,1)').^2+(N(:,2)-N(:,2)').^2;
m.a_ts=rand(n,n)*10;
m.a_to=m.a_ts+rand(n,n)*10;
m.v=rand(n,1)*40;
m.gamma_r=1;
m.gamma_p=1;
m.fixedprice=0.5;

[relocations,prices]=RelocationPricing(m)


