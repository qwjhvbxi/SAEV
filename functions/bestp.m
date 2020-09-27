
function p=bestp(dista,gamma_r,gamma_alt)

P=0:0.001:0.5;
R=@(p,d) (  exp(-d*p)./(exp(-d*p)+exp(-gamma_alt*d))  ).*(d*(p-gamma_r)); % net revenues at certain price and distance
[~,pind]=max(R(P,dista),[],2);
p=P(pind);

end