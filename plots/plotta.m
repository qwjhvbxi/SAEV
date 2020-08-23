function plotta(Res,PlotType,FolderName)

DataFolder=setDataFolder();
tsim=length(Res.Sim.relodist);
x=linspace(0,24,tsim);
x1=linspace(0,24,tsim+1);
x2=linspace(0,24,1440);
xt=0:4:24;
% priceUnits='(yen/kWh)';
priceUnits='(€/MWh)';
% Format='-depsc2';
% Resolution=[];
Format='-dpng';
Resolution='-r300';

if nargin<3
    FolderName=[];
end

switch PlotType
    
    case 'power'
        
        % power exchanged
        figure('Units','centimeters','Position',[10,7,10,7])
        plot(x,sum(Res.Sim.e,2)/1000)
        ylabel('power (MW)')
        hold on
        yyaxis right
        plot(x,Res.Params.elep(1:tsim),'r-')
        xlim([0,24])
        xticks(xt)
        xlabel('hours')
        ylabel(['electricity price ' priceUnits])
        legend({'power','price'},'Orientation','horizontal')
        set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
        if ~isempty(FolderName)
            print([DataFolder 'figures/' FolderName '/power'],Format,Resolution);
        end
        
end

