function prettyfigure(Half)

if nargin>0 && Half
    figure('Units','centimeters','Position',[10,7,7,5])
else
    figure('Units','centimeters','Position',[10,7,10,7])
end
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
box on

end