
%% check for constraints violations

convertvarformat

if ~exist('printplots','var')
    printplots=1;
end


if printplots==1
    
okok=1;
for t0=1:P.tsim %i
    
    eq3=round(sum(sum((p_check(:,:,:,t0))),3),5)<=1;     % (3)
    
    eq5=sum(u_check(:,:,t0))+sum(sum((p_check(:,:,:,t0))),3)==1;     % (5)
    
    eq6=sum(u_check(:,:,t0+1)+squeeze(sum(v_check(:,:,:,t0),2))+squeeze(sum(w_check(:,:,:,t0),2)))<=1;       % (6)
    
    eq8=round(sum(v_check(:,:,:,t0),3)-(d_check(:,:,t0)+c_check(:,:,t0)),5)<=0;   % (8)
    
    eq12=prod(  ones(P.n^2,1)*q_check(:,t0)'   >=  round(reshape(v_check(:,:,:,t0),[P.n^2,P.m])  *P.ad   .*(reshape(P.t,[P.n^2,1])*ones(1,P.m)) ,4) );
    
    eqNew=(abs(e_check(:,t0)) <= sum(u_check(:,:,t0+1))');
    
%     eqMacro=(squeeze(sum(sum(v_check(:,:,:,t0),1),2))  +  econ_check(:,t0) <=1);
    eqMacro=1; % no need to check this because it can vary with new alogrithm
    
    if(prod(eq3)*prod(eq5)*prod(eq6)*prod(eq8)*prod(eq12)*prod(eqNew)*prod(eqMacro)~=1)
        t0
        okok=0;
        if(prod(eq3)~=1)
            eq3
        end
        if(prod(eq5)~=1)
            eq5
        end
        if(prod(eq6)~=1)
            eq6
        end
        if(prod(eq8)~=1)
            eq8
        end
        if(prod(eq12)~=1)
            eq12
        end
        if(prod(eqNew)~=1)
            eqNew
        end
        if(prod(eqMacro)~=1)
            eqMacro
        end
    end
    
end
if okok==1
    fprintf('\n\nno constraints were violated! \n\n');
    clear eq*
end

end

if 0
    
   
    moving=squeeze(sum(sum(p_check(:,:,:,2:P.tsim+1),1),3));
    q_check(:,2:P.tsim+1)==q_check(:,1)+cumsum(e_check,2)-cumsum(e2_check,2)-moving*P.ad
    
    provarange=1:10;
    calcolato=q_check(:,provarange+1)
    derivato=q_check(:,1)+cumsum(e_check(:,provarange),2)-cumsum(e2_check(:,provarange),2)-moving(:,provarange)*P.ad
    
    calcolato-derivato
    
end



if 0
    
    movingvehicles=squeeze(sum(sum(p_check(:,:,:,:),1),3));
    idlevehicles=squeeze(sum(u_check(:,:,:),1));
    
    checkid=2;

    [car,distancetodest]=ind2sub([P.m,P.maxt+1],find(squeeze(sum(p_check(:,:,:,checkid),1))));
    disttodest=zeros(P.m,1);
    disttodest(car)=distancetodest;
            % soc               idle          moving        distance to destination 
    check1=[q_check(:,checkid) idlevehicles(:,checkid) movingvehicles(:,checkid) disttodest disttodest*P.ad]
    
    check2=[sum(chargeveh(1:P.m,:),2) sum(chargeveh(P.m+1:end,:),2)]
    
    [check1 check2 e_check(:,checkid-1) e2_check(:,checkid-1)]
    
    
    
    
    
    
    
    
%     q_check(:,checkid)+sum(chargeveh,2)-sum(chargeveh,2)
    
    fixedch
    chargeveh
    notmoving
    
    
    aeqplus
    
    
    beqplus
    




    d_checkm=reshape(X(1:P.n^2,:),[P.n,P.n,Pmacro.tsim+1]);

    p_checkm=reshape(X(P.n^2+1:P.n^2+P.n*P.m*(Pmacro.maxt+1),:),[P.n,P.m,Pmacro.maxt+1,Pmacro.tsim+1]);

    u_checkm=reshape(X(P.n^2+P.n*P.m*(Pmacro.maxt+1)+1:P.n^2+P.n*P.m*(Pmacro.maxt+2),:),[P.n,P.m,Pmacro.tsim+1]);

    q_checkm=reshape(X(P.n^2+P.n*P.m*(Pmacro.maxt+2)+1:P.n^2+P.n*P.m*(Pmacro.maxt+2)+P.m,:),[P.m,Pmacro.tsim+1]);

    c_checkm=reshape(Pmacro.c(1:P.n*P.n*Pmacro.tsim),[P.n,P.n,Pmacro.tsim]);

    v_checkm=reshape(zmacro(1:P.n*P.n*P.m,1:Pmacro.tsim),[P.n,P.n,P.m,Pmacro.tsim]);

    w_checkm=reshape(zmacro(P.n*P.n*P.m+1:2*P.n*P.n*P.m,1:Pmacro.tsim),[P.n,P.n,P.m,Pmacro.tsim]);

    e_checkm=zmacro(2*P.n*P.n*P.m+1:2*P.n*P.n*P.m+P.m,:);

    if P.v==3
        e2_checkm=zmacro(2*P.n*P.n*P.m+P.m+1:2*P.n*P.n*P.m+2*P.m,:);
    else
        e2_checkm=0;
    end

    charge_checkm=repelem(zmacro(end-P.m+1:end,1:Pmacro.tsim),1,P.macrosteps);

    econ_checkm=repelem(zmacro(end-P.m+1:end,1:Pmacro.tsim)~=0,1,P.macrosteps);


end




if 0

%% check objectives

namesim=['data/optimatrices-' num2str(P.thor) '-' num2str(P.n) '-' num2str(P.m) '-' num2str(P.v) '-' num2str(P.maxt)];
load(namesim);


fx;
fu; 
fq;

ttime=P.thor;
trialtime=500+(1:ttime);

% waiting
reshape(z(:,trialtime),size(z,1)*ttime,1)'*fx

% rebalancing
reshape(z(:,trialtime),size(z,1)*ttime,1)'*fu

% charging
reshape(z(:,trialtime),size(z,1)*ttime,1)'*(fq.*(repelem(P.elep(trialtime),P.ctrno+P.ctrnoplus,1)))

% final soc
reshape(z(:,trialtime),size(z,1)*ttime,1)'*fsoc

% charging fixed
% reshape(z(:,trialtime),size(z,1)*10,1)'*fq




namesimmacro=['data/optimatrices-' num2str(Pmacro.thor) '-' num2str(P.n) '-' num2str(P.m) '-' num2str(P.v) '-' num2str(Pmacro.maxt)];
load(namesimmacro);

ttime=Pmacro.thor;
trialtime=6+(1:ttime);

% waiting
waitcost=reshape(zmacro(:,trialtime),size(zmacro,1)*ttime,1)'*fx

% rebalancing
rebcost=reshape(zmacro(:,trialtime),size(zmacro,1)*ttime,1)'*fu

% charging
chargecost=reshape(zmacro(:,trialtime),size(zmacro,1)*ttime,1)'*(fq.*(repelem(P.elep(trialtime),P.ctrno+P.ctrnoplus,1)))

% discharging
dischcost=reshape(zmacro(:,trialtime),size(zmacro,1)*ttime,1)'*(fqv2g.*(repelem(P.elep(trialtime),P.ctrno+P.ctrnoplus,1)))

% final soc
soccost=reshape(zmacro(:,trialtime),size(zmacro,1)*ttime,1)'*fsoc

% charging fixed
% reshape(z(:,trialtime),size(z,1)*10,1)'*fq


waitcost*P.macrosteps
rebcost*Pmacro.rho1
-soccost*P.rho3
P.rho2*(chargecost+dischcost*0.9)


end


