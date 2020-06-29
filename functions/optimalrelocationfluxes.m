%% X=OPTIMALRELOCATIONFLUXES(F,R,T,limite)
% Caculate optimal relocation fluxes from imbalances
% 
% F [1 x N]: vector of excess vehicles at nodes
% R [1 x N]: vector of needed vehicles at nodes
% T [N x N]: matrix of travel times between nodes
% limite: time limit for relocation
% 
% Output X [N x N] is the matrix of optimal flows between nodes.
%
% See also FractionalRelocation2 stackablerelocation2

function X=optimalrelocationfluxes(F,R,T,limite)

% if there are imbalances and available vehicles
if sum(R)>0 && sum(F)>0
    
    if nargin<4
        limite=max(max(T)*2);
    end
    
    N=length(F); % stations
    
    % objective function (N^2 variables, all movements between stations)
    f=reshape(T,N*N,1)-repelem((R>0)',N,1)*limite;
    
    % constraints
    A1=[repmat(speye(N),1,N); kron(speye(N),ones(1,N))];
    b1=[F';R'];
    
    % integer variables (all)
    intcon=1:N*N;
    
    % all variables non-negative
    lb=zeros(N*N,1);
    
    try
        % launch optimization
        x=intlinprog(f,intcon,A1,b1,[],[],lb,[],optimoptions('intlinprog','display','none'));
    catch
        % octave formulation
        ctype=char('U'*ones(1,size(A1,1)));
        vtype=char('I'*ones(1,N*N));
        sense=1;
        [x,ov]=glpk(f,A1,b1,lb,[],ctype,vtype,sense);
    end
    
    % matrix form
    X=round(reshape(x,N,N));
    
else
    
    X=[];
    
end
    
end



