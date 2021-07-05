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

c=max(0.1,m.c);   % distance matrix (km)
c2=c(:);
n=size(c,1);      % number of nodes
prices=m.altp./c; % initial prices
K=m.K;            % fixed term (travel time cost)

% revenues=zeros(maxiter+1,1);
% ThisRevenue=CalcRevenue(m,prices);
% revenues(1)=ThisRevenue;


%% useful functions

% value at z
g=@(A,z) exp(-z.*c-K)./(exp(-z.*c-K)+A);

% derivative at z
d=@(A,z) -(A.*c.*exp(z.*c+K))./((1+A.*exp(z.*c+K)).^2); 


%% iterations

delta=zeros(n^2,m.maxiter);
priceshist=zeros(n^2,m.maxiter);

fprintf('\n iterations: ')

for k=1:m.maxiter

    % coefficients of y=Dx+C
    D=d(exp(-m.altp),prices);
    C=g(exp(-m.altp),prices)-D.*prices;
    p0=-C./D;
    
    % m.pmin=p1;
    % m.pmin=prices-(prices-p1)/k;
    m.pmin=prices-0.5./c; 

    if isfield(m,'fixedfleet') && m.fixedfleet
        m.pmax=p0;
    else
        m.pmax=prices+0.5./c;
    end
    
    m.amin=D.*m.pmax+C;%g(exp(-m.altp),-m.pmax.*c);
    m.amax=D.*m.pmin+C;%g(exp(-m.altp),-m.pmin.*c);
    
    Empty=(m.a==0);
    m.amin(Empty)=0;
    m.amax(Empty)=0;
    m.pmin(Empty)=0;
    m.pmax(Empty)=0;

    % launch pricing optimization for this iteration
    
    InputCheck=((m.pmax>=m.pmin).*(m.amax>=m.amin));
    if prod(InputCheck(:))==0
        warning('impossible')
        return
    end


    %% initializations

    a=m.a(:);
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
    % dreal=a(j)*exp(-p*c(j))./(exp(-p*c(j))+exp(-0.25*c(j)))

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
    
    % maxv=m.v;
    maxv=ones(n,1)*V;
    
    % bounds
    lb=[zeros(n^2,1);    pmin];
    ub=[repmat(maxv,n,1)*Relocation; pmax];
    ub(1:n+1:n^2)=0; % no relocation in same node

    if isfield(m,'fixedprice') && ~isempty(m.fixedprice)

        % bounds for fixed pricing
        lb(n^2+1:2*n^2)=m.fixedprice;
        ub(n^2+1:2*n^2)=m.fixedprice;

    end

    % cost function
    H=2*diag([zeros(n^2,1);a.*c2.*h]);
    f=[  c2*m.movingcostkm ;   -a.*c2.*(  s(:) + h*m.movingcostkm  )  ];

    options=optimoptions('quadprog','display','none');

    % optimization
    [X0,fval]=quadprog(H,f,A,b,[],[],lb,ub,[],options);

    % results (shift from convention [d,o] to [o,d])
    X1=round(X0,3);
    reloc=reshape(X1(1:n^2),n,n);
    newprices=reshape(X1(n^2+1:end),n,n);

    % there are no relocation actions
    if sum(reloc(:))==0
        reloc=[];
    end
    
    delta(:,k)=(newprices(:)-prices(:));
    priceshist(:,k)=newprices(:);
    
    prices=newprices;
    
    fprintf('*');
    
end

end


