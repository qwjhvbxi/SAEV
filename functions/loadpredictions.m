
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

switch double(P.Sim.mpcpredict)
    case 0 % legacy 
        k=strfind(P.tripfolder,'_');
        Pb.tripfolder=P.tripfolder(1:k(end)-1);
        Pb.tripday=P.tripday;
        Pbratio=str2double(P.tripfolder(k(end)+1:end)); % TODO: change!!
        [As2,Atimes2,~,~]=gettrips(Pb);
        [fo,fd]=computeod(As2,Atimes2,[],nc);
        Aforecast=matrixform(As2,Atimes2,nc,Pbratio);
    case 1 % perfect prediction
        [fo,fd]=computeod(As,Atimes,[],nc);
        Aforecast=matrixform(As,Atimes,nc,1);
    case 2 % aggregate+noise
        Aforecast=[];
        relativenoise=0.2;
        absolutenoise=0.1;
        [fo0,fd0]=computeod(As,Atimes,[],nc);
        fo=max(0,fo0.*(1+relativenoise*randn(size(fo0)))+absolutenoise*randn(size(fo0)));
        fd=max(0,fd0.*(1+relativenoise*randn(size(fd0)))+absolutenoise*randn(size(fo0)));
end

end

function Aforecast=matrixform(As,Atimes,nc,Pbratio)

% compute matrix form
maxt=max(Atimes(:));
As2ind=sub2ind([nc,nc],As(:,1),As(:,2));
Aforecast=sparse(Atimes(:,1),As2ind,1,maxt,nc^2)/Pbratio;

end




