%% LANCIAGA2(P,K,O)
% launch genetic optimization for parameters tx, ts, tr.
% P is the struct with the basic parameters. K is the fleet size, O is the
% number of operators. If given as ranges (two values), the GA finds the
% best fleet size and no. of operators. In this case, economic information
% are needed [NOT YET IMPLEMENTED FOR RANGES].

function lanciaga2(P,K)

if nargin==2
    % fixed K
    if length(K)==1
        K(2)=K;
    end
else
    % K part of optimization (need to specify other objective fun.)
    warning('need to specify different objective (economic)')
    K=[5000 15000];
end

save('tempP.mat','P');

%% variable's boundaries
%   tx ts tr K
lb=[1  5  5  K(1)];
ub=[30 30 30 K(2)];
intcon=1:4;

% options=optimoptions('ga','PlotFcn', @gaplotbestf);
options=gaoptimset('PopulationSize',60,'Display','iter','UseParallel','Always','OutputFcns',@dispfunc);

[x,fval,exitflag,output]=ga(@fitnessga,5,[],[],[],[],lb,ub,[],intcon,options);

save('resultGA2.mat');

end

