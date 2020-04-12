%% 

[~,R1]=generateplotline2(1,'dropped','O',Ovec,'reloctype',2,'K',5000,'scenarioid',1:10);
[~,R2]=generateplotline2(1,'dropped','O',Ovec,'reloctype',2,'K',5000,'eta',2,'scenarioid',1:10);
[~,R3]=generateplotline2(1,'dropped','O',0,   'reloctype',4,'K',5000,'scenarioid',1:10);
[~,R4]=generateplotline2(1,'dropped','O',0,   'reloctype',3,'K',5000,'scenarioid',1:10);
[~,R5]=generateplotline2(1,'dropped','O',OvecShort2,'reloctype',7,'K',5000,'scenarioid',[7 9 10],'t',60,'ts',60,'tr',60);
% [~,R5]=generateplotline2(1,'dropped','O',OvecShort2,'reloctype',7,'K',5000,'scenarioid',9,'t',60,'ts',60,'tr',60);
% [~,R6]=generateplotline2(1,'dropped','O',OvecShort,'reloctype',6,'K',5000,'scenarioid',1:10,'Ib',100); % Boyaci
[~,R6]=generateplotline2(1,'dropped','O',round(Ovec/2),'reloctype',6,'K',5000,'scenarioid',1:10,'Ib',100); % Boyaci
[~,R7]=generateplotline2(1,'dropped','O',OvecShort,'reloctype',2,'K',5000,'eta',1,'scenarioid',1:10);

firstel=1;
firstel2=2;
lw=1;
c=lines(7);
figure('pos',[500 200 400 300])
hold on
[p1,~]=errbarmedian(Ovec(firstel:end),R1(firstel:end,:)'*100,lw,1,c(1,:)); % OpR-E8
[p2,~]=errbarmedian(Ovec(firstel:end),R2(firstel:end,:)'*100,lw,1,c(2,:)); % OpR-E3
[p3a,p3b]=errbarmedian(Ovec(firstel:end),R3*100,lw,1,c(3,:)); % UsR
[p4a,p4b]=errbarmedian(Ovec(firstel:end),R4*100,lw,1,c(4,:)); % AR
R5b=repelem(R5,1,[5 3 2]);
[p5,~]=errbarmedian(OvecShort2(firstel:end),R5b(firstel:end,:)'*100,lw,1,c(5,:)); % TrR real
% p5=plot(OvecShort2(firstel:end),R5(firstel:end)*100,'LineWidth',lw,'Color',c(5,:),'Marker','x');  % TrR
% [p5,~]=errbarmedian(OvecShort2(firstel2:end),([R5(firstel2:end,:)-0.02 R5(firstel2:end,:) R5(firstel2:end,:)+0.02]')*100,lw,0,c(5,:)); % TrR
[p6,~]=errbarmedian(Ovec(firstel:end),R6(firstel:end,:)'*100,lw,1,c(6,:)); % Boyaci
% [p7,~]=errbarmedian(OvecShort(firstel:end),R7(firstel:end,:)'*100,lw,1,c(7,:)); % OpR-E2
[lgd, icons, plots, txt]=legend([p1 p2  p3b p4b p5 p6],cases{[1:2,4:end]});
xlabel('relocators')
ylabel('dropped requests (%)')
