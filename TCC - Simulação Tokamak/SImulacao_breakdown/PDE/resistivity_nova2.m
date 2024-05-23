function out = resistivity_nova2(ne,Te,Zeff)

e = 1.6022e-19; me = 9.1094e-31; eo = 8.8542e-12;

Fz = (1 + 1.198*Zeff + 0.222*Zeff.^2)./(1 + 2.966*Zeff + 0.753*Zeff.^2);
Ln = log(12*pi*eo^1.5/e^1.5*Te.^1.5./sqrt(ne));

out.eta_par = sqrt(2*me*e)/(12*pi^1.5*eo^2)*Zeff.*Fz.*Ln./Te.^1.5;
out.eta_perp = 1.96*out.eta_par;
out.v_ei = e^2/me*ne.*out.eta_par;

end