% 
% INPUT PARAMETERS
% 
%           seedno: '20190225150309'
%      arrivalseed: 132
%      eleprofseed: 4
%      eleproftype: 'T'
%       mpcpredict: 0
%             type: 1
%       macrosteps: 15
%             thor: 12
%                e: 2
%      tripsperday: 600
%             rho1: 0.01
%             rho2: 0.0001
%             rho3: 0.002
%             rho4: 1e-06
%                n: 10
%                m: 20
%          battery: 50
%            delta: 1.5
%         speedkmh: 20
%         chargekw: 20
%      consumption: 0.15
%      cyclingcost: 10
%       initialsoc: 0.7
%           minsoc: 0.2
%           maxsoc: 1
%        v2gminsoc: 0.5
%           mapmod: 'Small'
% 
% 
% WORKING PARAMETERS
%                 
%             tsim: 720
%            speed: 0.66667
%               ac: 0.013
%               ad: 0.002
%                v: 2
%         stations: [10×2 double]
%     tripsperhour: [24×1 double]
%             elep: [1440×1 double]
%                c: [73200×1 double]
%        cexpected: [73200×1 double]
%                t: [10×10 double]
%             maxt: 12
%            varno: 2900
%        varnoplus: 20
%            ctrno: 4000
%        ctrnoplus: 20
%                x: [2920×1 double]



Pe.vehiclecost=5e6/5/365;   % yen/day/vehicle
Pe.batterycost=20000;       % yen/kWh
Pe.chargercost=5e5/10/365;  % yen/day/charger (level II)
Pe.gasolineprice=150;       % yen/liter
Pe.ICconsumption=5/100;     % liter/km
Pe.cyclelife=2000;          % battery cycle life

% total cost of charging/moving
Y.totalcost=(P.type>0)*P.battery*(sum(e_check*P.elep(1:P.tsim))-0.9*sum(e2_check*P.elep(1:P.tsim))  ) ...
           +(P.type==0)*Y.totalkm*Pe.ICconsumption*Pe.gasolineprice;
Y.totalcostfinalsoc=(P.type>0)*P.battery*(sum(e_check*P.elep(1:P.tsim))-0.9*sum(e2_check*P.elep(1:P.tsim)) -( sum(q_check(:,end))-sum(q_check(:,1)) )*mean(P.elep(1:P.tsim))  ) ... % mean
           +(P.type==0)*Y.totalkm*Pe.ICconsumption*Pe.gasolineprice;

% total fixed costs + battery costs
Y.fixedcost=Pe.vehiclecost+  (P.type>0)*(Pe.batterycost/Pe.cyclelife*(sum(e_check*P.elep(1:P.tsim))- (sum(q_check(:,end))-sum(q_check(:,1)) ) ) ...
                             +sum(max(squeeze(sum(u_check,2)),[],2))*Pe.chargercost);







