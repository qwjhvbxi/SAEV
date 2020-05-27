

addpath plots

[P2,R2]=generateplotline3('NYC2016',[],'Operations.maxwait',20,'modechoice',[false true],'m',[8000 9000 10000]);

[~,R2d]=generateplotline3('NYC2016','dropped','Operations.maxwait',20,'modechoice',[false true],'m',[8000 9000 10000]);


for j=1:3
    Walking(j)=1-sum(R2(2,j).Sim.chosenmode)/length(R2(2,j).Sim.chosenmode);
    Served=((R2(2,j).Sim.chosenmode==1).*(R2(2,j).Sim.dropped==0));
    Profits(j,2)=
end

figure('Units','centimeters','Position',[10,7,10,7])
hold on
plot(R2d'*100,'o-')
ylabel('unserved %')
yyaxis right
plot(Walking*100,'s-')
ylabel('%')
legend({'no mc';'mc';'walking'})




%% waiting time estimation vs real

R1=R2(2,3);
Selection=logical((R1.Sim.chosenmode==1).*(R1.Sim.dropped==0).*(R1.Sim.waitingestimated>0).*(R1.Sim.waiting>0));
prova=full([R1.Sim.waitingestimated(Selection) R1.Sim.waiting(Selection)]);
M=zeros(10,10);
for j=1:10
    for k=1:10
        M(j,k)=sum(prod(prova==[j k]*2,2));
    end
end
figure
imagesc(2:2:20,2:2:20,M)

