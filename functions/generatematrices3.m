% [nomefile]=GENERATEMATRICES3(n,m,Tr,maxt,ac,ad,thor,maxsoc,minsoc)
% generate optimization matrices 
% 
% see also GENERAL11

function [nomefile]=generatematrices3(n,m,Tr,maxt,ac,ad,thor,maxsoc,minsoc)

DataFolder=setDataFolder();
hashcell={n,m,Tr,maxt,ac,ad,thor,maxsoc,minsoc,'v3'};
nomefile=[DataFolder 'optimatrices/' DataHash(hashcell) '.mat'];

if ~exist(nomefile,'file')

%% preliminaries

v=3;
varno=n^2+m*n*(2+maxt); % passengers waiting, cars positions, cars waiting |  ~58 million in my model
varnoplus=(v>0)*m;   % soc
ctrno=n^2*m*2;
ctrnoplus=(v>=1)*m+(v==3)*m;


%% vector d,p,u,q

dveclength=n*n;%length(d_vec);
pveclength=n*m*(maxt+1);%length(p_vec);
uveclength=n*m;%length(u_vec);
qveclength=m*(v>=1);

chargecontrolslength=m*((v>=1)+(v==3));
blength=n*n*m+n*n*m+chargecontrolslength;
xlength=dveclength+pveclength+uveclength+qveclength;

t_vec=reshape(Tr,[n*n,1]);


%% creation of coefficient matrices (eq. 10) [d p u q]

fprintf('creation of coefficient matrices \n');


% d         d2=A1*Param.x+B1*b+c2; % eq.1

A1=[diag(ones(n^2,1)) zeros(n^2,n*m*(maxt+2))];
B1=[repmat(diag(-ones(n^2,1)),1,m)  ,   zeros(n^2 ,    n^2*m   +  chargecontrolslength)  ];

fprintf(' -> d done \n');


% p         p2=A2*Param.x+B2*b; % eq 2

A2=[false(n*m*(maxt+1),n^2) diag(true(n*m*(maxt),1),n*m) false(n*m*(maxt+1),n*m)];
B2=zeros(n*m*(maxt+1),n*n*m);
for jtime=1:maxt+1
    for jcar=1:m
        for jarrival=1:n
            % select all the columns where arrival station =j and time to
            % station from departure to arrival is =z-1+1
            %           previous cars                                       current car                                         next cars
            B2((jtime-1)*n*m+(jcar-1)*n+jarrival,:)   =   [zeros(1,n^2*(jcar-1)) ,  zeros(1,n*(jarrival-1)) ,  (Tr(jarrival,:)==(jtime)) ,   zeros(1,n*(n-jarrival)) ,  zeros(1,n^2*(m-jcar))];
        end
    end
end
B2=[B2 B2 zeros(size(B2,1),chargecontrolslength)];

fprintf(' -> p done \n');


% u          u2=A3*Param.x+B3*b; % eq 4

A3=zeros(n*m,n*m*(maxt+1));
A3(1:n*m+1:(n*m)^2)=1;          % p part
A3=[zeros(n*m,n^2) A3 diag(ones(n*m,1))];
B3=[];
for j1=1:m
    B3temp=[];
    for j2=1:n
        B3temp=[B3temp diag(ones(n,1))];
    end
    B3=[B3;zeros(n,n^2*(j1-1)) , B3temp , zeros(n,n^2*(m-j1))];
end
B3=-[B3 B3 zeros(size(B3,1),chargecontrolslength)];

fprintf(' -> u done \n');


% q

B4=[];
A4=[];
if v>=1
    
    %           q % eq new
    selectorp=logical(repmat(kron(eye(m),ones(n,1)),maxt+1,1));
    
    % p
    for j1=1:m
        A4=[A4; -ad*sum(A2(selectorp(:,j1),:))];
        B4=[B4; -ad*sum(B2(selectorp(:,j1),:))];
    end
    
    % q
    A4=[A4 eye(m)];
    
    % v,w,e
    B4(:,n^2*m*2+1:end)=B4(:,n^2*m*2+1:end)  +   [eye(m) , -eye(m*(v==3))]; % check!!!
    
    fprintf(' -> q done \n');
    
end


Aopti=sparse([[[sparse(A1);sparse(A2);sparse(A3)] sparse(varno,m*(v>=2))];sparse(A4)]);
Bopti=sparse([B1;B2;B3;B4]);
clear A1 A2 A3 B1 B2 B3

fprintf(' -> all done \n\n');

% mattests


%% generation of constraints matrices for optimization

fprintf('creation of matrices for optimization problem \n');

% objective function for x
fx=zeros(blength*thor,1);
fsoc=zeros(blength*thor,1);

% objective function for secondary objective (minimize redistribution of vehicles)
fu=repmat(  [   zeros(n^2*m,1) ; repmat(reshape(Tr,[n^2,1]),[m,1])   ; zeros(chargecontrolslength,1) ]   ,   [thor,1]  );

% objective function for charging
fq=zeros(blength*thor,1);
fqv2g=0;%zeros(blength*thor,1);
if v>=1
    fq=repmat([zeros(n^2*m*2,1) ; ones(m,1) ; zeros(m*(v==3),1)],thor,1);
    if v==3
        fqv2g=repmat([zeros(n^2*m*2,1) ; zeros(m,1) ; -ones(m,1)],thor,1);
    end
end

% number of constraints needed
constq=(v>0)*m   +  (v>=1)*m*2 + (v==3)*m*2;  % number of constraints for charging (e1, e2) need to be connected to charging station

% generation of matrices for (dis)equality constraints
Aeq=sparse([]);
beq=sparse([]);
beqc=sparse([]);
Adis=[];
bdis=[];
bdisc=[];
bdisC=[];

fprintf(char('_'*ones(1,thor)))
fprintf('\n');

for j=1:thor        % j = tau
    
    
    %% generate x(t+1) for current thorizon
    
    AB0=[];
    for s=1:j
        AB0=[AB0 Aopti^(j-s)*Bopti];
    end
    AB0=[AB0 sparse(size(Bopti,1),size(Bopti,2)*(thor-j))];
    Ac0=Aopti^j;
    AC0=sparse([[repmat(eye(n^2,'logical'),1,j) sparse(n^2,n^2*(thor-j))];sparse(size(Bopti,1)-dveclength,thor*n^2)]);    % to multiply with C vector size(Bopti,1) x n^2*thor
    
    
    
    
    %% === minimization function === %%
    
    % generate vector fx (as in fx*u) (main objective)
    
    fx=fx+sum(AB0(1:dveclength,:))';
    
    
    % generate vector fsoc (as in fsoc*x) (secondary objective)
    
    fsoc=fsoc+sum(AB0(varno+1:varno+varnoplus,:))';
    
    
    %% === constraints === %%
    
    %% generate Aopti
    
    %% constraints on x (7) (0 or 1 / or from 0 to 1)
    
    Adis=[Adis; sparse( AB0(dveclength+1:end,:) )   ; sparse(  -AB0(1:end,:)) ];
    bdis=[bdis; sparse([   [   -Ac0(dveclength+1:end,:)  ;   Ac0(1:end,:) ]  , sparse(2*(varno+varnoplus)-dveclength,varnoplus*(v<1)) ]) ];
    bdisc=[bdisc; sparse([ ones(varno-dveclength,1) ; ones(varnoplus,1)*maxsoc     ;   sparse(varno,1) ; -ones(varnoplus,1)*minsoc])  ;  sparse(constq,1) ];
    bdisC=[bdisC; sparse([   -AC0(dveclength+1:end,:)  ;   AC0(1:end,:) ])  ;  sparse(constq,size(AC0,2)) ];
    
    Adistemp1=[];
    bdistemp1=[];

    for jcar=1:m

        %% constraints for enough charge for a trip

        Adistemp1=[Adistemp1;  [   zeros(1,(j-1)*blength) ...        % previous timesteps
            repmat([sparse(1,n^2*(jcar-1))  sparse(t_vec'*ad)+minsoc  sparse(1,n^2*(m-jcar))],[1,2])  ... % v, w (energy needed for each trip)
            zeros(1,chargecontrolslength) ...
            zeros(1,(thor-j)*blength) ] ... %sparse(1,(thor-j-1)*blength)  ] ...    % next timesteps
            +  (    -AB0(varno+jcar,:)      )];    % q
        bdistemp1=[bdistemp1; (   Ac0(varno+jcar,:) )];  % q


        %% constraints on connection to charging station

        % first constraint (positive values)
        Adistemp0=sparse(-sum(AB0(dveclength+pveclength+(jcar-1)*n+1:dveclength+pveclength+jcar*n,    :))); % u
        Adistemp2=Adistemp0;
        Adistemp2(1,(ctrno+ctrnoplus)*(j-1)+ctrno+jcar)=Adistemp2(1,(ctrno+ctrnoplus)*(j-1)+ctrno+jcar)+1;  % e
        Adistemp1=[Adistemp1;Adistemp2];
        bdistemp1=[bdistemp1; sum(Ac0(dveclength+pveclength+(jcar-1)*n+1:dveclength+pveclength+jcar*n,:)) ]; % u

        % second constraint (negative values)
        Adistemp2=Adistemp0;
        Adistemp2(1,(ctrno+ctrnoplus)*(j-1)+ctrno+jcar)=Adistemp2(1,(ctrno+ctrnoplus)*(j-1)+ctrno+jcar)-1;  % e
        Adistemp1=[Adistemp1;Adistemp2];
        bdistemp1=[bdistemp1; sum(Ac0(dveclength+pveclength+(jcar-1)*n+1:dveclength+pveclength+jcar*n,:))]; % u


        %% constraints on connection to charging station for V2G
        
        if v==3 

            % first constraint (positive values)
            Adistemp2=Adistemp0;
            Adistemp2(1,(ctrno+ctrnoplus)*(j-1)+ctrno+m+jcar) =   Adistemp2(1,(ctrno+ctrnoplus)*(j-1)+ctrno+m+jcar)+1;  % e2
            Adistemp1=[Adistemp1;Adistemp2];
            bdistemp1=[bdistemp1; sum(Ac0(dveclength+pveclength+(jcar-1)*n+1:dveclength+pveclength+jcar*n,:))]; % u

            % second constraint (negative values)
            Adistemp2=Adistemp0;
            Adistemp2(1,(ctrno+ctrnoplus)*(j-1)+ctrno+m+jcar) =   Adistemp2(1,(ctrno+ctrnoplus)*(j-1)+ctrno+m+jcar)-1;  % e2
            Adistemp1=[Adistemp1;Adistemp2];
            bdistemp1=[bdistemp1; sum(Ac0(dveclength+pveclength+(jcar-1)*n+1:dveclength+pveclength+jcar*n,:))]; % u

        end

    end

    Adis=[Adis;Adistemp1];
    bdis=[bdis;bdistemp1];
        
    
    %% generate Aeq
    
    % constraints on vii wii
    
    tmp2=sparse(n,n*n);
    tmp2(1:n*n+n+1:n^3)=1;
    tmp30=sparse(n*m,n*n*m);
    tmp3=tmp30;
    
    for jk=1:m
        
        tmp3(n*(jk-1)+1:n*jk,n*n*(jk-1)+1:n*n*jk)=tmp2;
        
    end
    
    tmp4=[tmp3 ,  tmp30 , sparse(size(tmp3,1),chargecontrolslength)  ; tmp30 ,  tmp3  , sparse(size(tmp3,1),chargecontrolslength)];
    
    tmp5=[sparse(size(tmp4,1),(j-1)*blength)  , tmp4 ,  sparse(size(tmp4,1),(thor-j)*blength) ];
    
    Aeq=[Aeq;tmp5];
    
    beq=[beq;zeros(size(tmp5,1),size(Aopti,2))];
    
    beqc=[beqc;zeros(size(tmp5,1),1)];
    
    
    fprintf(char('*'))
    
end


%% generate ub lb

lb=zeros(n*n*m*2*thor,1);

ub=ones(n*n*m*2*thor,1);

if (v>=1)
    lb=repmat( [zeros(n*n*m*2,1) ;  zeros(chargecontrolslength,1) ] , thor,1); % lower bound changed !!!
    ub=repmat( [ones(n*n*m*2,1)  ;  ones(chargecontrolslength,1)*ac ] , thor,1);
end


%% intcon

intcon=[];
for j3=1:thor
    intcon=[intcon (j3-1)*(ctrno+ctrnoplus)+(1:ctrno)];
end


%% save matrices

fprintf('\n\n');

save(nomefile,'Aopti','Bopti','Adis','bdis','bdisc','bdisC','Aeq','beq','beqc','intcon','lb','ub','fu','fx','fq','fqv2g','fsoc','-v7.3');

end

end









