%% [u,r,d,e,C]=scenarioreduction(x[,p,C,u])
% Calculate new reduced scenarios u with associated equivalent probability
% r from original scenarios x with original probability p. 
% C is the matrix of distances between scenarios.
% Specify p, C, u if known
% Leave p=[] or do not give input if all are equal probability.
% Each column of x is a scenario. 
% d is the matrix containing the index of the new scenario closest to the
% the original scenario.
% e is the matrix containing the distance between the closest new scenario
% to the the original scenario. 

function [u,r,d,e,C]=scenarioreduction(x,p,C,u)
    
    %% calculations
    
    n=size(x,2); % number of original scenarios
    L=size(x,1); % cardinality of scenarios
    
    if nargin==1 || isempty(p)
        p=ones(n,1)/n;  % probability of scenarios
    end
    
    if nargin==1
        
        %% calculate distances between scenarios (C)

        addpath('../emd')

        % normalize profiles
        pnorm=(x-min(min(x)))/max(max(sum(x-min(min(x)))));

        % calculate distances with EMD
        C=zeros(n);         % distance matrix
        D=ones(L)-diag(ones(L,1));  % auxiliary matrix for calculations
        for i=1:n
            for j=1:n
                C(i,j)=emd_mex(pnorm(:,i)'/sum(pnorm(:,i)),pnorm(:,j)'/sum(pnorm(:,j)),D);
            end

            if rem(i,10)==0
                clc
                n-i
            end
        end
    
    end
    
    
    if nargin<3
    
        %% forward selection algorithm (u)

        s=n;

        u=zeros(s,1);   % selected scenarios
        q=zeros(s,1);   % new probability of selected scenarios

        z=p'.*sum(C,'omitnan');
        [~,u(1)]=min(z,[],'omitnan');

        for i=2:s

            % create function z'
            for oi=1:n

                % remove omega' (oi) from results
                c2=C;
                c2(oi,:)=NaN;
                c2(:,oi)=NaN;

                z(oi)=p(oi)*sum(  min(c2( : , u(u~=0) ),[],2,'omitnan') ,'omitnan' );

            end

            z(u((u~=0)))=NaN; % do not consider already selected solutions

            [~,u(i)]=min(z,[],'omitnan');

            if rem(i,10)==0
                clc
                n-i
            end

        end
    
    end
    
    
    %% calculate new probabilities
    
    r=zeros(n,n);
    d=zeros(n,n);
    e=zeros(n,n);
    
    for i=1:n
        
        [r(:,i),d(:,i),e(:,i)]=newprobs(C,p,u,i);
        
        if rem(i,10)==0
            clc
            n-i
        end
        
    end
    
    clc
    fprintf('done.\n\n')

end




    
    






%% plots

function plotta

    % plot image of resulting matrix
    figure
    imagesc(C)

    % plot most different/similar pairs of scenarios
    dso=reshape(tril(C),n*n,1);
    [~,Is]=sort(dso,'descend');
    for i=[1, 2, round(n^2/2-n/2)]
        [a,b]=ind2sub([n,n],Is(i));
        figure
        plot([pday(:,a) pday(:,b)])
    end

end










% %% first algorithm
% s=5;
% c=(Distances);
% u=zeros(s,1);
% p=ones(365,1)/365;
% J=zeros(365,1);
%
% % z=p'.*sum(c);
% % [~,u(1)]=min(z);
% % J(u(1))=[];
% % [J;u1];
% % c(u(1),:)=NaN;
% % c(:,u(1))=NaN;
%
% for i=1:s
%
%     z=p'.*sum(c,'omitnan');
%     z(u((u~=0)))=NaN;
%
%     [~,u(i)]=min(z,[],'omitnan');
%
%
%     J(u(i))=NaN;
%     c(u(i),:)=NaN;
%     c(:,u(i))=NaN;
%
%     % per ogni k, prendi il minimo tra uk e u1k
%     c=min(c,c(:,u(i))*ones(1,365));
%
% end


%
% z=p'.*sum(c,'omitnan');
% [~,u(1)]=min(z,[],'omitnan');
%
%
% for i=2:s
%
%
% %     J(u(i))=NaN;
% %     c(u(i),:)=NaN;
% %     c(:,u(i))=NaN;
%
%     % per ogni k, prendi il minimo tra uk e u1k
% %     c=min(c,c(:,u(i))*ones(1,365));
%
%     for oi=1:365 % j = o'
%         c2=c;
%         c2(oi,:)=NaN;
%         c2(:,oi)=NaN;
%
%         z(oi)=sum(  min(c2( : , u(u~=0) ),[],2,'omitnan')  );
%
%     end
% %     [~,u(i)]=min(res)
%
%
%     z(u((u~=0)))=NaN;
%
%     [~,u(i)]=min(z,[],'omitnan');
%
% end



% %% second algorithm (sbagliato)
% s=5;
% c=(Distances);
% u=zeros(s,1);
% p=ones(365,1)/365;
% n=365;
% i=1;
% while n>s && i<=length(Iq)
%
%     [a,b]=ind2sub([365,365],Iq(i));
% %     figure
% %     hold on
% %     plot([pday(:,a) pday(:,b)])
%     % prendo solo a e salvo, aggiungo la probabilita di b
%     p(a)=p(a)+p(b);
%     p(b)=0;
%
%     n=sum((p~=0));
%     i=i+1;
%
% end
% u=find(p);



