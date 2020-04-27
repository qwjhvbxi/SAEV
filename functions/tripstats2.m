



function [Atimes,fo,fd,dk]=tripstats2(A,Atimes,T)

N=size(T,1);

% cleanup A
Atimes(Atimes(:,1)==0,1)=1;
Atimes(Atimes(:,2)==0,2)=1;

fo=zeros(double(max(Atimes(:,2)))*2,N);
fd=zeros(double(max(Atimes(:,2)))*2,N);
dk=zeros(double(max(Atimes(:,2)))*2,1);

for i=1:1440
    ThisMinute=logical(Atimes(:,1)==i);
    fo(i,:)=accumarray([A(ThisMinute,1);N],[ones(sum(ThisMinute),1);0]);
    fd(i,:)=accumarray([A(ThisMinute,2);N],[ones(sum(ThisMinute),1);0]);
    dk(i)=sum(T(sub2ind(size(T),A(ThisMinute,1),A(ThisMinute,2))));
end

fo=sparse(fo);
fd=sparse(fd);

end