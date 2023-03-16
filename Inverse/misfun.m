function[f,df,d2f] = misfun(Rho,Rho_R)
%[f,df,d2f] = misfun(t)
%
t   = Rho - Rho_R;
f   = 0.5*(t'*t);
df  = t;
d2f = ones(numel(t),1);

return

