%% C=fitga(x)
% calculate fitness (cost) for a specific individual

function C=fitga(x)

load('tempGA.mat','P')
P.charging=repmat((x>0)'.*x',2,1);
P.discharging=repmat((x<0)'.*x',2,1);
R=general11(P,-1,0);
C=R.cost+sum(R.Sim.waiting(:,3))*10;

end