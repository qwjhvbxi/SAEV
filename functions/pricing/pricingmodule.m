function [perDistanceTariff,surcharges,altp]=pricingmodule(Pricing,forecastOD,alternativeCosts,ui)

n=Pricing.n;
    
% initialize matrix of fare per minute
perDistanceTariff=ones(n,n).*Pricing.basetariff;

% initialize surcharges per stations
surcharges=zeros(2*n,1);

if ~isempty(forecastOD)

    % if Pricing.dynamic>0 && dynamicpricing<3 && (i==1 || mod(i-(ts+tr+1),tp)==0)
    if Pricing.dynamic 

        nodebased=isfield(Pricing,'nodebased') && Pricing.nodebased;

        % number of vehicles at each station (including vehicles directed there)
        Pricing.v=histc(ui,1:n);
        a_tp=sparse(forecastOD(:,1),forecastOD(:,2),1,n,n);%+q_t;
        a_tp(1:n+1:end)=0;
        Pricing.a=a_tp;

        if numel(Pricing.alternative)>1
            [a,Ib,~]=unique(forecastOD,'rows','stable');
            altp=sparse(a(:,1),a(:,2),alternativeCosts(Ib),n,n);
        else
            altp=Pricing.alternative.*Pricing.c; % alternative price for each OD
        end
        Pricing.altp=altp;

        if ~nodebased

            [perDistanceTariff,~,~]=NLPricing5(Pricing); % OD-pair-based pricing

        else 

            [surcharges,~,~]=NLPricingNodes(Pricing); % node-based pricing

        end

    else

        % expected trips
        if numel(Pricing.alternative)>1
            [a,Ib,~]=unique(forecastOD,'rows','stable');
            altp=sparse(a(:,1),a(:,2),alternativeCosts(Ib),n,n);
        else
            altp=Pricing.alternative.*Pricing.c; % alternative price for each OD
        end

    end

else
    
    altp=zeros(n,n);
    
end

end