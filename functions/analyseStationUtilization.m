addpath ../model_aev5/functions/




T=P.T;

% Duration=15;

Resolution=10;
Inizio=12*30/Resolution;

RealVehiclesAtNode=histc(Res.Sim.u(1:Resolution:end,:)',1:n);
% RealVehiclesAtNode=histc(Res.Sim.u',1:n);
VehiclesAtNode=RealVehiclesAtNode(:,1:Inizio);
EffectiveVehiclesAtNode=VehiclesAtNode;
TripMatrix=squeeze(sum(reshape(Res.Trips.c,n,n,Resolution,length(Res.Trips.c)/Resolution),3));
% prova
% TripMatrix=[0 3 1;1 0 2;3 4 0];
% T=rand(3,3);
% T=tril(T)+tril(T)';
% T(1:size(T,1)+1:end)=0;
% n=size(TripMatrix,1);
% TripMatrix(1:n+1:end)=0;



Origins=squeeze(sum(TripMatrix,2));
Destinations=squeeze(sum(TripMatrix));

for i=Inizio:1440/P.e/Resolution

% MinutesSelection=(i-1)*30+1:(i-1)*30+Duration/P.e;
% TripMatrix=sum(Res.Trips.c(:,:,MinutesSelection),3);

% TripMatrix=Res.Trips.c;
% MinutesIdle=histc(Res.Sim.u',1:n);

% Origins(:,i)=squeeze(sum(TripMatrix,2));
% Destinations(:,i)=squeeze(sum(TripMatrix)');

ExtraVehicles=VehiclesAtNode(:,i)+Destinations(:,i)-Origins(:,i);
MinExtra=0;

% identify optimal relocation flux
F=ExtraVehicles.*(ExtraVehicles>MinExtra);
R=-(ExtraVehicles-MinExtra).*(ExtraVehicles<MinExtra);
x=optimalrelocationfluxes(F',R',T,100);

RelocationOrigins(:,i)=sum(x,2);
RelocationDestinations(:,i)=sum(x)';

VehiclesAtNode(:,i+1)=VehiclesAtNode(:,i)+Destinations(:,i)-Origins(:,i)+RelocationDestinations(:,i)-RelocationOrigins(:,i);
EffectiveVehiclesAtNode(:,i+1)=VehiclesAtNode(:,i+1)-(Destinations(:,i)+Origins(:,i)+RelocationDestinations(:,i)+RelocationOrigins(:,i))/10;


end

return

k=40
figure
hold on
VError=VehiclesAtNode(:,k)-RealVehiclesAtNode(:,k);
plot(VError)
% plot(Destinations(:,k))
% plot(Origins(:,k))

scatter(VError,Destinations(:,k-1))%+Origins(:,k-1))

figure
hold on
plot(sum((RealVehiclesAtNode-VehiclesAtNode)),'r')
plot(sum(abs(RealVehiclesAtNode-VehiclesAtNode)),'r--')
plot(sum((RealVehiclesAtNode-EffectiveVehiclesAtNode)),'b')
plot(sum(abs(RealVehiclesAtNode-EffectiveVehiclesAtNode)),'b--')




figure
imagesc(VehiclesAtNode)
figure
imagesc(RealVehiclesAtNode)

figure
hold on
for k=1:24
    clf
    plot(VehiclesAtNode(:,k))
    hold on
    plot(RealVehiclesAtNode(:,k))
    drawnow
    % scatter(VehiclesAtNode(:,k),RealVehiclesAtNode(:,k))
    % line([0 300],[0 300])
    pause(0.5)
end


% [e,F]=emd_mex(Destinations',Origins',T);
% F=[F(:,1:2)+1 F(:,3)];
% F((F(:,1)==F(:,2)),:)=[];
% F=[F;[n n 0]];

% RelocationOrigins=accumarray(F((F(:,1)~=0),1),F((F(:,1)~=0),3));
% RelocationDestinations=accumarray(F((F(:,2)~=0),2),F((F(:,2)~=0),3));

% Origins2=Origins+RelocationOrigins;
% Destins2=Destinations+RelocationDestinations;
% 
% % Origins2=Origins+RelocationDestinations;
% % Destins2=Destinations'+RelocationOrigins;
% 
% [Destinations Origins RelocationOrigins RelocationDestinations]
% 
% [Origins2 Destins2]




%% 

MinutesIdle=histc(reshape(Res.Sim.u(MinutesSelection,:),length(MinutesSelection)*P.m,1),1:n);
VehiclesAtNode=histc(Res.Sim.u(i*30,:),1:n)';

[VehiclesAtNode Origins2 Destins2 MinutesIdle (VehiclesAtNode+Destins2./Origins2)*Duration/P.e]


figure
plot([MinutesIdle VehiclesAtNode*Duration/P.e  (VehiclesAtNode+Destins2./Origins2)*Duration/P.e])

figure
plot([MinutesIdle (VehiclesAtNode)*Duration/P.e])


[VehiclesAtNode Origins2 Destins2 VehiclesAtNode+Destins2-Origins2 histc(Res.Sim.u(i*30+30,:),1:n)']

figure
hold on
plot(VehiclesAtNode,'--')
plot(VehiclesAtNode+Destins2-Origins2)
plot(histc(Res.Sim.u(i*30+30,:),1:n)')






figure
plot(Origins)
hold on
plot(Destinations)
yyaxis right
plot(MinutesIdle)

figure
scatter(MinutesIdle,Origins)

figure
scatter(MinutesIdle,Destinations)


figure
scatter(MinutesIdle,Destinations-Origins)




