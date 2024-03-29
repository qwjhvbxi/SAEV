%% [distancebasedtariff,surcharges,prices]=PRICINGMODULE(Pricing,forecastOD,alternativeCosts,ui)
% Launch pricing module optimization.
%
% See also: mainsim

function [distancebasedtariff,surcharges,prices]=pricingmodule(Pricing,forecastOD,alternativeCosts,ui)

% TODO: output should just be OD pricing, either calculated directly for
% each OD, or by distance+inboud/outboud node

n=size(Pricing.c,1);
    
% initialize matrix of fare per minute
distancebasedtariff=ones(n,n).*Pricing.basetariffkm;

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
        Pricing.altp=alternativeCosts;
        Pricing.K=Pricing.T*Pricing.VOT/60*Pricing.traveltimecost;

        if ~nodebased

            [distancebasedtariff,~,~]=nlpricingod(Pricing); % OD-pair-based pricing

        else 

            [surcharges,~,~]=nlpricingnode(Pricing); % node-based pricing

        end

    end

end

surchargemat=surcharges(1:n)+surcharges(n+1:2*n)';
prices=max(Pricing.mintariff,distancebasedtariff.*Pricing.c)+surchargemat;

end