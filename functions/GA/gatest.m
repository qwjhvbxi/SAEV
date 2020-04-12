
params
P.enlayeralg='aggregate';
P.trlayeralg='simplified';
R=general11(P,-1,1)




P.enlayeralg='custom';
save('tempGA.mat','P')
% P.charging=rand(48*2,1);
% P.discharging=zeros(48*2,1);
% R=general11(P,-1,-1)

figure
hold on
N=48;
lb=-1*ones(N,1);
ub=ones(N,1);
options=gaoptimset('PopulationSize',200,'Display','iter','UseParallel','Always','OutputFcns',@dispfunc);%,'PlotFcn','gaplotbestf');%
[x,fval,exitflag,output]=ga(@fitga,N,[],[],[],[],lb,ub,[],[],options);

