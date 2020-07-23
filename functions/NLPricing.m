


function [prices,k]=NLPricing(m)

f=@(a,s) exp(s)./(exp(s)+a);
d=@(a,s) a*exp(s)./((a+exp(s)).^2);  % derivative at s

% linearization intervals (7 intervals, 8 limits)
Points=f(exp(0),-3.5:3.5);

n=size(m.c,1);  % nodes
p2=-m.gamma_p;   % base price
a=exp(p2*m.c); % exp of benefit of alternative modes for each node pair

m.w=zeros(n,n);    % position of linearization (between -3 and 3)
u=0.5;
maxIter=4;

for k=1:maxIter

    % limits to normalized price (corresponding to between 0 and 2*baseprice)
    m.pmin=(p2*m.c+m.w+u)./(2*p2*m.c);
    m.pmax=(p2*m.c+m.w-u)./(2*p2*m.c);
    m.pmax(1:n+1:end)=1;
    m.pmin(1:n+1:end)=0;

    % probabilities of trip acceptance at the limits
    m.amin=1-Points(5+m.w);
    m.amax=1-Points(4+m.w);

    [~,prices]=RelocationPricing3(m);

    movedown=(prices==m.pmin);
    moveup=(prices==m.pmax);

    if sum(movedown(:))+sum(moveup(:))>0
        m.w=m.w-movedown;
        m.w=m.w+moveup;
    else
        break;
    end

end


% fmin=f(a,2*p2*pmin0.*m.c);
% fmax=f(a,2*p2*pmax0.*m.c);
% 
% fmin=f(a,p2*m.c+w-u);
% fmax=f(a,p2*m.c+w+u);
% pmax=1-fmin;
% pmin=1-fmax;
% pmax(1:n+1:end)=0.5;
% pmin(1:n+1:end)=0.5;
% 
% m.pmin=
% m.pmax=


return


%% testing 

N=round([   0.65574 , 0.70605
      0.03571 , 0.03183
      0.84913 , 0.27692
      0.93399 , 0.04617
      ]*30);
A=[ 0 3 2 0;
    1 0 9 8;
    1 1 0 0;
    3 0 6 0];
m.a=A;
% m.v=[5;18;2;9];
m.v=[5;0;2;9];

m.c=(N(:,1)-N(:,1)').^2+(N(:,2)-N(:,2)').^2;
m.gamma_r=0.1;
m.gamma_p=0.5;
% m.fixedprice=0.5;

prices=NLPricing(m)
