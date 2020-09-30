%% [X,prices]=RelocationPricing5(m)
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

function [X,prices]=RelocationPricing6(m)

InputCheck=((m.pmax>=m.pmin).*(m.amax>=m.amin));
if prod(InputCheck(:))==0
    warning('impossible')
    return
end

n=size(m.c,1);    % nodes

if isfield(m,'relocation')
    Relocation=double(logical(m.relocation));
else
    Relocation=1;
end
m.amin(1:n+1:end)=0;
m.amax(1:n+1:end)=1;
m.pmin(1:n+1:end)=-10;
m.pmax(1:n+1:end)=10;


%% initializations

q=m.a;
q(1:n+1:end)=0;

d=m.c;
Tariff=m.gamma_d(d+eye(n));
Da=m.amax-m.amin;
Dp=m.pmax-m.pmin;


%% demand function
% demand is in the form d(:)=a*(s-h*p)

alpha=m.amax+(Da./Dp).*m.pmin; % static term of demand
beta=Da./Dp; % p multiplier for demand


%% constraints

% relocation [n^2]   demand [n^2]   surcharge,discount o [n*2]      surcharge,discount d [n*2]

% constraint on relocation
Ar=repmat(eye(n),1,n)-kron(eye(n),ones(1,n));
A=[Ar , Ar , sparse(n,n*4)];
b=m.v;

% constraints on relocation vehicles
Av=[repmat(eye(n),1,n) sparse(n,n^2+n*4)];
A=[A;Av];
b=[b;m.v];

% constraints for linearization
Al=[sparse(n^2,2*n^2) , repmat(eye(n),n,1) , -repmat(eye(n),n,1) , repelem(eye(n),n,1) , -repelem(eye(n),n,1) ];
A=[A;Al;-Al];
b=[b;m.pmax(:);-m.pmin(:)];

% constraint on demand
Aeq=[];
beq=[];
M=q(:).*beta(:);
A=[A;[sparse(n^2,n^2) , speye(n^2) , repmat(eye(n),n,1).*M , -repmat(eye(n),n,1).*M , repelem(eye(n),n,1).*M , -repelem(eye(n),n,1).*M ]];
b=[b;q(:).*alpha(:)];

% bounds
lb=zeros(2*n^2+n*4,1);
ub=[repmat(m.v*Relocation,n,1); q(:); ones(n*4,1)*5];
%ub=[repmat(m.v,n,1); q(:); ones(n*2,1)*5 ; zeros(n*2,1)];
ub(1:n+1:n^2)=0; % no relocation in same node

if isfield(m,'fixedprice') && ~isempty(m.fixedprice)

    % bounds for fixed pricing
    lb(n^2+1:2*n^2)=m.fixedprice;
    ub(n^2+1:2*n^2)=m.fixedprice;
    
end

% cost function
M=1;
% H0=spdiags([zeros(n^2,1);sum(beta.*q,2)+0.001;(sum(beta.*q))'+0.001],0,2*n^2+n*4,2*n^2+n*4);
H0=sparse(2*n^2+n*4,2*n^2+n*4);
H0(n^2+1:2*n^2,2*n^2+1:end)=[  -repmat(speye(n),n,1).*M , repmat(speye(n),n,1).*M , -repelem(speye(n),n,1).*M , repelem(speye(n),n,1).*M  ];

% H0(2*n^2+1:2*n^2+n,2*n^2+n+1:2*n^2+2*n)=speye(n);
% H0(2*n^2+2*n+1:2*n^2+3*n , 2*n^2+3*n+1:2*n^2+4*n)=speye(n);
% H0(2*n^2+1:2*n^2+n,2*n^2+2*n+1:2*n^2+3*n)=speye(n);
% H0(2*n^2+n+1:2*n^2+2*n , 2*n^2+3*n+1:2*n^2+4*n)=speye(n);

H=(H0+H0');
f=[ d(:)*m.gamma_r ; ... relocation
    -d(:).*(Tariff(:)-m.gamma_r); ... demand
    ones(n*4,1)*0.0001];

if 0
    options=optimoptions('fmincon','display','iter','MaxFunctionEvaluations',5e4); 
    myobjfun = @(x)parameterfun(x,H,f);
    x0=[zeros(1,2*n^2),max(0,m.pvec(1:n)),-min(0,m.pvec(1:n)),max(0,m.pvec(n+1:2*n)),-min(0,m.pvec(n+1:2*n))];
    [X0,fval,Exitflag] = fmincon(myobjfun, x0, A, b, Aeq, beq, lb, ub,[],options);%, @myg)
else
    options=optimoptions('quadprog','display','none');
    % optimization
    [X0,fval,Exitflag]=quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
end


if Exitflag>=0

    X1=round(X0(:),3);
    X=reshape(X1(1:n^2),n,n);
    pricevec=X1(2*n^2+1:end);
    prices=[pricevec(1:n)-pricevec(n+1:2*n) ; pricevec(2*n+1:3*n)-pricevec(3*n+1:4*n)];
    demand=reshape(X1(n^2+1:2*n^2),n,n);
    % there are no relocation actions
    if sum(X(:))==0
        X=[];
    end

else
    
    X=[];
    prices=[];
    
end

return


%% debug

% m.v=[3 0 0]';
% m.c=[0 2 3
%      2 0 4
%      3 4 2];
% m.a=[0 10 4
%      3  0 9
%      7  1 0];

m.v=[0 0 20]';
m.c=[0 2 3
     2 0 4
     3 4 0];
m.a=[0 10 0
     5  0 2
     0  0 0];
%  m.relocation=0;

% m.pmin=rand(3)/2;
% m.pmax=rand(3)/2+0.5; 
% m.amin=rand(3)/2;
% m.amax=rand(3)/2+0.5;
m.pmin=ones(3)*0.25;
m.pmax=ones(3)*0.75;
m.amin=ones(3)*0.25;
m.amax=ones(3)*0.75;
m.gamma_r=0.1;

[reloc,prices,demand]=RelocationPricing4(m)



