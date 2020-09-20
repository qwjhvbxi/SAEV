%% [X,prices]=RelocationPricing4(m)
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

function [X,prices]=RelocationPricing4(m)

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


%% constraints

% constraint on relocation
Ar=repmat(eye(n),1,n)-kron(eye(n),ones(1,n));
A=[Ar , sparse(n,n^2) , Ar];
b=m.v;

% constraint on demand
Ad=[sparse(n^2,n^2)  ,  spdiags(a.*h,0,n^2,n^2)  ,  speye(n^2)];
A=[A;Ad];
b=[b;a.*s(:)];

% constraints on relocation vehicles
Av=[repmat(eye(n),1,n) sparse(n,n^2*2)];
A=[A;Av];
b=[b;m.v];

% bounds
%   relocations      prices   demand
lb=[zeros(n^2,1);    pmin;    zeros(n^2,1)];
ub=[repmat(m.v,n,1); pmax;    inf(n^2,1)];
ub(1:n+1:n^2)=0; % no relocation in same node

if isfield(m,'fixedprice') && ~isempty(m.fixedprice)

    % bounds for fixed pricing
    lb(n^2+1:2*n^2)=m.fixedprice;
    ub(n^2+1:2*n^2)=m.fixedprice;
    
end

HShape=[0 0 0;0 0 1;0 1 0];
H=kron(HShape,spdiags(-c,0,n^2,n^2) );
f=[  c*m.gamma_r ; zeros(n^2,1) ;  c.*m.gamma_r  ];

options=optimoptions('quadprog','display','none');

% optimization
[X0,fval]=quadprog(H,f,A,b,[],[],lb,ub,[],options);

% results (shift from convention [d,o] to [o,d])
X1=round(X0,3);
X=reshape(X1(1:n^2),n,n);
prices=reshape(X1(n^2+1:2*n^2),n,n);
demand=reshape(X1(2*n^2+1:3*n^2),n,n);

% there are no relocation actions
if sum(X(:))==0
    X=[];
end

return


% %% debug
% 
% m.v=[3 0 0]';
% m.c=[0 2 3
%      2 0 4
%      3 4 2];
% m.a=[0 10 4
%      3  0 9
%      7  1 0];
% m.pmin=rand(3)/2;
% m.pmax=rand(3)/2+0.5;
% m.amin=rand(3)/2;
% m.amax=rand(3)/2+0.5;
% m.gamma_r=0.1;
% 
% RelocationPricing4(m)



