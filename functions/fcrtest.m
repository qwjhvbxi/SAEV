function [Res]=fcrtest(P,R,plotta)

if isfield(P,'FCR') && ~isempty(P.FCR)

    FileName=P.FCR.filename;
    GridDay=P.gridday;
    Limits=P.FCR.limits;
    Contracted=P.FCR.contracted;
    Power=sum(R.Sim.ef,2);

    DataFolder=getdatafolder();
    [fraw,~,fresolution]=readexternalfile([DataFolder 'grid/' FileName '.csv'],GridDay,false);
    F=average2(fraw,P.Sim.e/fresolution);

    FCRe=min(1,max(-1,(F-50)/(Limits(2)-Limits(1))*2)); % needed FCR (EV POV)
    FCRe=FCRe(1:1440);
    PowerExpect=FCRe*Contracted*1000;

    Res.DeltaPower=sparse(double(round(Power-PowerExpect)));
    Res.FailMinutes=(sum(abs(Res.DeltaPower)>0));
    Res.Totals=[mean(max(0,Power))*24 , mean(max(0,-Power))*24];

    if nargin>2 && plotta
        figure
        hold on 
        plot(PowerExpect)
        plot(Power)
    end

else
    
    Res.DeltaPower=nan;
    Res.FailMinutes=nan;
    Res.Totals=nan;
    
end

end