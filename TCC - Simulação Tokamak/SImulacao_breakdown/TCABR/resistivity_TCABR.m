function out = resistivity_TCABR(ne,Te,Zeff)
global G e me
if Te==0 Te=G.Te0; end 
if ne == 0 ne=G.n0; end
eo = 8.8542e-12;

Fz = (1 + 1.198*Zeff + 0.222*Zeff.^2)./(1 + 2.966*Zeff + 0.753*Zeff.^2);
Ln = log(12*pi*eo^1.5/e^1.5*Te.^1.5./sqrt(ne));

out.eta_par = sqrt(2*me*e)/(12*pi^1.5*eo^2)*Zeff.*Fz.*Ln./Te.^1.5;
out.eta_perp = 1.96*out.eta_par;
out.v_ei = e^2/me*ne.*out.eta_par;

end