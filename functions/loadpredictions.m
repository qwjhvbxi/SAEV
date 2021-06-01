
% TODO: cleanup call to secondary trip file (real vs expected/forecasted)
% for simulations with modechoice, OD prediction is needed, but can be
% probability distributions. For other simulations, aggregates are enough
% 2 options: 1) OD probability distributions in matrix form
%            2) aggregates at nodes in the form fo,fd
% option 2) can be easily calculated from 1) and used by relocation,
% pricing etc.
% option 1) only for pricing

function [fo,fd,Aforecast]=loadpredictions(P,As,Atimes)

nc=max(As(:));

if ~P.Sim.mpcpredict
    k=strfind(P.tripfolder,'_');
    Pb.tripfolder=P.tripfolder(1:k(end)-1);
    Pb.tripday=P.tripday;
    Pbratio=str2double(P.tripfolder(k(end)+1:end)); % TODO: change!!
    [As2,Atimes2,~,~]=gettrips(Pb);
%     cumulativeTripArrivals2=cumulativeTripArrivals2(1:P.Sim.e:end);
else
    Pbratio=1;
    As2=As;
    Atimes2=Atimes;
%     cumulativeTripArrivals2=cumulativeTripArrivals;
end

% compute matrix form
maxt=max(Atimes(:));
As2ind=sub2ind([nc,nc],As(:,1),As(:,2));
Aforecast=sparse(Atimes,As2ind,1,maxt,nc^2)/Pbratio;

% compute fo and fd
[fo,fd]=computeod(As2,Atimes2,[],nc);



