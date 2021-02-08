
function plotmap(A,c)
s=scatter(A(:,1),A(:,2),1,'filled'); 
s.MarkerFaceColor=c;
s.MarkerFaceAlpha=0.05;
axis equal tight
box on
% set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
end