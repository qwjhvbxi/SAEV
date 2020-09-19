
function p=bestp(c)

P=0:0.001:0.5;
R=@(p,c) (  exp(-p.*c)./(exp(-p.*c)+exp(-0.25*c))  ).*(p-0.1).*c; % net revenues at certain price and distance
[~,pind]=max(R(P,c));
p=P(pind);

end