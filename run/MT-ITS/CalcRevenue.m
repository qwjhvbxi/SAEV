
function ThisRevenue=CalcRevenue(m,prices,reloc)

p1=prices.*m.c;
Aeff=( exp(-p1)./(exp(-p1)+exp(-m.gamma_alt*m.c))  ).*m.a;
if ~isempty(reloc)
    relocation=sum(sum(reloc.*m.c))*m.gamma_r;
    r=sum(reloc)'-sum(reloc,2);
else 
    relocation=0;
    r=0;
end
ThisRevenue=full(sum(sum(Aeff.*(p1-m.gamma_r*m.c))))-relocation;

s=m.v+sum(Aeff)'-sum(Aeff,2);
% s+r