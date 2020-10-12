function plotta(Res,PlotType,FolderName,FileID,Format,Resolution)

DataFolder=setDataFolder();
tsim=length(Res.Sim.relodist);
x=linspace(0,24,tsim);
x1=linspace(0,24,tsim+1);
x2=linspace(0,24,1440);
xt=0:4:24;
% priceUnits='(yen/kWh)';
% priceUnits='($/MWh)';
priceUnits=['(' char(8364) '/MWh)']; % €
if nargin<6
    if nargin<5
        % Format='-depsc2';
        % Resolution=[];
        Format='-dpng';
        Resolution='-r300';
    else 
        Resolution=[];
    end
end

if nargin<3
    FolderName=[];
end

if nargin<4
    FileID='';
end

switch PlotType
    
    case 'power'
        
        if isfield(Res.Sim,'ef')
            ef=double(sum(Res.Sim.ef,2));
        else
            ef=zeros(size(Res.Sim.e,1),1);
        end
        
        % power exchanged
        figure('Units','centimeters','Position',[10,7,10,7])
        stairs(x,full(sum(Res.Sim.e,2)+ef)/1000)
        ylabel('power (MW)')
        hold on
        yyaxis right
        stairs(x,Res.Params.elep(1:tsim),'r-')
        xlim([0,24])
        xticks(xt)
        xlabel('hours')
        ylabel(['electricity price ' priceUnits])
        legend({'power','price'},'Orientation','horizontal')
        set(gca,'FontUnits','points','FontWeight','normal','FontName','Times')
        if ~isempty(FolderName)
            print([DataFolder 'figures/' FolderName '/power' FileID],Format,Resolution);
        end
        
end

