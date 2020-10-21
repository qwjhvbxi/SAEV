
DataFolder=setDataFolder();

Pav=repmat([4 3 4.5 -4 -3.5 -4.8],2,1);
Pset=repmat([4.7 0 -3],2,1);
Pop=repmat([4 0 -3],2,1);
Pfcr=repmat(repelem(Pop(1,:),1,10)+randn(1,30)/5,2,1);
x=[0:10:20;10:10:30];
x2=[0:29;1:30];
figure('Units','centimeters','Position',[10,7,10,4])
hold on
C=lines(10);
line(x,Pop,'LineStyle','-','Color',C(3,:),'LineWidth',2);
line(repmat(x,1,2),Pav,'LineStyle','--','Color',C(1,:).^(1/3),'LineWidth',2);
line(x,Pset,'LineStyle',':','Color',C(2,:),'LineWidth',2);

line(x2(:,1:26),Pfcr(:,1:26),'LineStyle','-','Color',C(3,:),'LineWidth',2);

line([26,26;27,27]',[-5,5;-5,5]','LineStyle',':','Color',[0 0 0],'LineWidth',1.5)
xlabel('time steps')
grid on
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
% print([DataFolder 'figures/FAU/diagram1'],'-depsc2');

%%

Pav=repelem([4 3 4.5;-4 -3.5 -4.8],1,10)';
Pset=repelem([4.7 0 -3],1,10)';
Pop=repelem([4 0 -3],1,10)';
Pfcr=Pop+randn(30,1)/5;

DataFolder=setDataFolder();
figure('Units','centimeters','Position',[10,7,10,7])
hold on
line([1,length(Pav)],[0 0],'LineStyle',':');

stairs(Pop,'-','LineWidth',1.5)
stairs(Pset,'-.','LineWidth',1)
stairs(Pav,'k--','LineWidth',1.5)
stairs(Pfcr,'-','LineWidth',1)
xlabel('time steps')
grid on
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
% print([DataFolder 'figures/FAU/diagram1'],'-depsc2');

