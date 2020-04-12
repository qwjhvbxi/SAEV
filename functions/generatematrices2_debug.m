%% debug
addpath('functions')
addpath('utilities')
n=10;
m=20;
v=3;
thor=12;
Tr=randi(11,n,n);%[0 3 2;4 0 1;3 3 0];
maxt=max(max(Tr));
ac=0.1;
ad=0.2;
maxsoc=1;
minsoc=0.3;
varno0=n^2+m*n*(2+maxt); % passengers waiting, cars positions, cars waiting |  ~58 million in my model
varnoplus=(v>0)*m;   % soc
ctrno=n^2*m*2;
ctrnoplus=(v>=1)*m+(v==3)*m;

tic
[nomefile]=generatematrices2(n,m,Tr,maxt,ac,ad,thor,maxsoc,minsoc);
toc
% 
% Ptemp.n=n;
% Ptemp.m=m;
% Ptemp.v=v;
% Ptemp.t=Tr;
% Ptemp.maxt=max(max(Ptemp.Tr));
% Ptemp.ac=ac;
% Ptemp.ad=ad;
% Ptemp.thor=thor;
% Ptemp.maxsoc=maxsoc;
% Ptemp.minsoc=minsoc;
% Ptemp.varno=Ptemp.n^2+Ptemp.m*Ptemp.n*(2+Ptemp.maxt); % passengers waiting, cars positions, cars waiting |  ~58 million in my model
% Ptemp.varnoplus=(Ptemp.v>0)*Ptemp.m;   % soc
% Ptemp.ctrno=Ptemp.n^2*Ptemp.m*2;
% Ptemp.ctrnoplus=(Ptemp.v>=1)*Ptemp.m+(Ptemp.v==3)*Ptemp.m;
% 
% tic
% generatematrices
% toc


