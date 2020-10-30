%
%
%
%
% test 1 instance of the problem.

%% generate scenario

DataFolder=setDataFolder();
addpath functions/pricing
addpath run/MT-ITS
P=cpar('NYC2018');
[A,Atimes,AbuckC,Distances]=loadTrips(P);
load([DataFolder 'scenarios/' P.scenario '.mat'],'T','Clusters','chargingStations');

h=8;
trips=AbuckC((h-1)*60+1):AbuckC((h-1)*60+30);
a=A(trips,:);
n=max(A(:));
a_tp=sparse(A(trips,1),A(trips,2),1,n,n);
a_tp(1:n+1:end)=0;

%% parameters

N=3;
V=1000; % number of vehicles
Iter=5;

Par.gamma_r=0.1; % relocation cost per minute
Par.gamma_p=0.25; % base tariff per minute
Par.gamma_alt=0.25;
Par.VOT=15;
Par.pricingwaiting=0;
Par.c=T;
Par.relocation=1;
% Par.v=zeros(n,1);
Par.a=a_tp;
Par.maxIter=Iter;

revenues2=zeros(N,Iter);

for k=1:N
    
    ui=randi(n,V,1);
    Par.v=histc(ui,1:n);
    
    %% without pricing optimization

    % prices=ones(n,n)*Par.gamma_p;
    Par1=Par;
    Par1.pmin=zeros(n);
    Par1.pmax=ones(n);
    Par1.amin=0.5;
    Par1.amax=0.5;
    Par1.fixedprice=Par1.gamma_p;

    [reloc,prices1,fval]=RelocationPricing7(Par1);

    revenues0=CalcRevenue(Par1,prices1,reloc);
    revenues1(k)=revenues0;


    %% with pricing optimization

    Par2=Par;
    [prices2,~,~,~,revenues]=NLPricing5(Par2); % OD-pair-based pricing
    
    revenues2(k,:)=revenues;


end

%% comparison


figure('Units','centimeters','Position',[10,7,10,7])
hold on
errbarmedian((1:Iter),revenues1',1,1,'b')
errbarmedian((1:Iter),revenues2,1,1,'k')
% errbarmedian((1:Iter),revenues3,1,1)
xlim([1,Iter]);
xticks(1:Iter)
xlabel('iterations')
ylabel('net revenues')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');

return

%% comparison and zoom

figure('Units','centimeters','Position',[10,7,10,7])
hold on
subplot(2,1,1)
[p1,p2]=errbarmedian((1:Iter),revenues1',1,1,'b')
errbarmedian((1:Iter),revenues2,1,1,'k')
box on
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
% xticklabels({})
xticks([])
ylabel('net revenues')
subplot(2,1,2)
errbarmedian((1:Iter),revenues2,1,1,'k')
% errbarmedian((1:Iter),revenues3,1,1)
xlim([1,Iter]);
xticks(1:Iter)
xlabel('iterations')
set(gca,'FontUnits','points','FontWeight','normal','FontSize',11,'FontName','Times');
print([DataFolder 'figures/Pricing/iterations'],'-depsc2');

%% plot example

figure
hold on

prices=[0.25,0.4,0.5,0.55];
% prices=[0.25,0.2,0.15,0.12];
% prices=[0.25,0.1,0.15,0.12];
m.gamma_alt=0.25;
m.c=10;
g=@(a,s) exp(s)./(exp(s)+a);        % value at s
d=@(a,z,c) -(a.*c.*exp(z.*c))./((1+a.*exp(z.*c)).^2); % derivative at s

for k=1:length(prices)

    % coefficients of y=Dx+C
    D=d(exp(-m.gamma_alt.*m.c),prices(k),m.c);
    C=g(exp(-m.gamma_alt.*m.c),-prices(k).*m.c)-D.*prices(k);
    
    % alternative
    p1=(1-C)./D;
    p0=-C./D;
    
%     m.pmin=p1;
    m.pmin=prices(k)-(prices(k)-p1)/k; % should be used
    m.pmax=p0;
%     m.pmax=prices+(p0-prices)/k; 
    
    m.amin=D.*m.pmax+C;%g(exp(-m.gamma_alt.*m.c),-m.pmax.*m.c);
    m.amax=D.*m.pmin+C;%g(exp(-m.gamma_alt.*m.c),-m.pmin.*m.c);
    
%     Empty=(m.a==0);
%     m.amin(Empty)=0;
%     m.amax(Empty)=0;
%     m.pmin(Empty)=0;
%     m.pmax(Empty)=0;

    j=1;
    line([m.pmin(j),m.pmax(j)],[m.amax(j),m.amin(j)]);
    p=0:0.01:1;
    plot(p,g(exp(-m.gamma_alt*m.c(j)),-p*m.c(j)),'k:')
    plot(prices(k),D*prices(k)+C,'o')
    plot(prices(k+1),D*prices(k+1)+C,'x')
    
end



