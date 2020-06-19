
function RelocationPricing()

% k=10;   % fleet size
% pmax=1; % max price
% c=[0 2 3;2 0 2;3 2 0]/10; % distance matrix
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
c=(N(:,1)-N(:,1)').^2+(N(:,2)-N(:,2)').^2;

n=size(c,1);    % nodes



% a_ts=[  0 0 0;
%         0 0 0;
%         0 0 0]; % demand up to ts
% a_to=[  0 6 0;
%         0 0 0;
%         9 0 0]; % demand up to to (should be always higher than previous)
% v=[0;10;0]; % vehicles at nodes

% a_ts=zeros(n,n);
% a_to=zeros(n,n);

a_ts=rand(n,n)*10;
a_to=a_ts+rand(n,n)*10;
v=rand(n,1)*10;


% variables: r, p 

% constraint on relocation
a_ts_v=reshape(a_ts,n^2,1);
a_to_v=reshape(a_to,n^2,1);
a0to=kron(eye(n),ones(1,n));
a0ts=repmat(eye(n),1,n);
a2=zeros(size(a0to));
a2(logical(a0to))=-a_to_v;
a2(logical(a0ts))=a2(logical(a0ts))+a_ts_v;
A=[kron(eye(n),ones(1,n))-repmat(eye(n),1,n)  ,  a2];
b=v+sum(a_ts,2)-sum(a_to,1)';

% % constraint on relocation backup
% a0=kron(eye(n),ones(1,n));
% a2=zeros(size(a0));
% a2(logical(a0))=a;
% A=[kron(eye(n),ones(1,n))-repmat(eye(n),1,n)  ,  a2];
% b=sum(a,1)'-v;

% bounds
lb=zeros(n^2*2,1);
ub=[repelem(v,n,1);ones(n^2,1)];
ub(1:n+1:n^2)=0;

H=2*diag([zeros(n^2,1);a_to_v]);

f=[reshape(c,n^2,1);-a_to_v];


X=quadprog(H,f,A,b,[],[],lb,ub);

X1=round(X,5);
relocations=reshape(X1(1:n^2),n,n)
prices=reshape(X1(n^2+1:end),n,n)
a_ts.*(1-prices)
a_to.*(1-prices)


% P=0:0.1:1;
% alf=10;
% figure
% hold on
% plot(P,P*alf-P.^2*alf)



end
