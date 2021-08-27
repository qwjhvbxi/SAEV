function [uinew,selected]=movetocs(Par,ui,gi,qi)

uinew=ui;
ncs=length(Par.cspos);

% initial status
mat1=(ui==Par.cspos);
atChargingStation=sum(mat1);
whichcs=(1:ncs)*mat1;

IdleReached=find((gi>0).*(1-atChargingStation));

% number of empty spots at each charging station
R=Par.cssize-histc(whichcs,1:ncs)';
fullCS=find(R<0);
R=max(0,R);
for k=1:length(fullCS)
    thisv=find(whichcs==fullCS(k));
    [~,y]=sort(qi(thisv));
    IdleReached=[IdleReached,thisv(y(Par.cssize(fullCS(k))+1:end))];
end

V=IdleReached;

% variables xij vehicle i to station j
% sum_i xij <= cssize_j
X=Par.Tr(uinew(V),Par.cspos)+(1-qi(V)');

nvar=length(V)*ncs;
intcon=1:nvar;
f=X(:)-30;
A=[repelem(speye(ncs),1,length(V)) ; repmat(speye(length(V)),1,ncs)];
b=[R;ones(length(V),1)];

lb=zeros(nvar,1);
ub=ones(nvar,1);

options=optimoptions('intlinprog','display','none');

x0=intlinprog(f,intcon,A,b,[],[],lb,ub,[],options);
x1=reshape(x0,length(V),ncs);
x=x1*Par.cspos;
selected=V(logical(x));
uinew(selected)=x(logical(x));

end
