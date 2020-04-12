

%% load/create reduced scenarios

% scenarioname='DayaheadPrices2017Germany';
% load(['data/' scenarioname])
% phour=DayaheadPrices2017Germany;
% x=reshape(phour,24,365);

scenarioname='DayaheadPrices2016TokyoNew';
load('../data/JEPX');
TokyoDAFY2016=(Jepx.M((365*3+1)*48+1:end,8));
x=reshape(TokyoDAFY2016,48,365);

try
    
    load(['../data/' scenarioname '-r']);

catch
    
    [u,r,d,e,C]=scenarioreduction(x);%,C,u);
    
    save(['../data/' scenarioname '-r'],'u','r','d','e','C','x');

end


%% analysis of options

% decreasing total error with more scenarios
figure
plot(sum(e));

% if I select only the first s scenarios from i reducted
s=[5 10 20];
Sh=zeros(n,length(s));
Er=zeros(n,length(s));

for k=1:length(s)
    for i=s(k):n

        % probability of reduced scenario
        q=r(1:i,i);

        [~,b]=sort(q,'descend');
        cm=cumsum(q(b));
        Sh(i,k)=cm(s(k)); % share of scenario represented (others are considered outliers)

        for j=1:s(k)
            Er(i,k)=Er(i)+sum(  e( d(:,i)==b(j) , i )   );
        end

    end
end

figure
hold on
for i=1:length(s)
    plot(s(i):n,1-Sh(s(i):n,i))
    text(100+s(i),1-Sh(100+s(i),i),num2str(s(i))) 
end
ylabel('excluded scenarios')
yyaxis right
plot(s(1):n,Er(s(1):n))
ylabel('error')
xlabel('no. reduced scenarios')
plot(sum(e),'--');


%% save selected scenarios

s=10;
r(1:s,s)
figure
hold on
for i=1:s

    profday=x(:,u(i));
    elep=repmat(profday,2,1); 
%     elep=repmat(repelem(profday,2,1),2,1); % for hourly profiles

    %     plot(pday(:,u(i)))%,'LineWidth',1+log(dayselection(u(i))))
%     plot(profday)

    save(['../data/elep/elep-T' num2str(u(i))],'elep')

end


