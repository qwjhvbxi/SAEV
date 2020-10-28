% function for GA
function y = parameterfun(x,H,f)
y=.5*x*H*x' + x*f;