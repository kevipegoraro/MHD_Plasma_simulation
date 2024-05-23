function breakdown

%% Input
Ntime = 1000;
dt = 1e-7;
R0 = 0.62;
B0 = 1.07;
Vloop = 7;
p0 = 0.005;
n0 = 1e4;
Te0 = 0.026;
gas = 'H2';

% Loading Green's functions table and plasma poloidal flux distribution
load('/home/kevi/greens_table/greens_table/green_table_65x65.mat');

% Initializing physical constants and matrices
e = 1.6e-19;
me = 9.11e-31;
mi = 1.67e-27;
[zg,rg] = meshgrid(z,r); dr = mean(diff(r)); dz = mean(diff(z));
Br   = zeros(size(rg)); Bz = Br;
n   = n0*ones(size(rg,1),size(rg,2),Ntime+1); pe = e*n*Te0; pion = pe;
Jr = ones(size(rg,1),size(rg,2),Ntime+1); Jz = Jr; Jphi = Jr;

% Current moments
Iq   = [1 1 0 0 0 0 -1 -1 -0.6 -0.6 0.6 0.6 0]*1e4;  % Quadrupole field

% Computing quadrupole and toroidal magnetic and electric fields
for icoil = 1:length(Iq)
    eval(['Br = Br + F' num2str(icoil) '.BR*Iq(' num2str(icoil) ');'])
    eval(['Bz = Bz + F' num2str(icoil) '.BZ*Iq(' num2str(icoil) ');'])
end
Bphi = R0*B0./rg;
Ephi = Vloop/2/pi./rg;

% Computing quadrupole and toroidal magnetic and electric fields
out = resistivity(n0,Te0,1);
a = first_townsend_coeff(p0,Ephi,gas,0);
Leff = R0*B0./(0.001+(Br.^2 + Bz.^2).^0.5);
rho = out.eta_par;
v_en = 7.89e11*p0;
v_ion = 0;
v_loss = 0;
v_ei = out.v_ei;

for it = 1:Ntime
    %keyboard
    % Electron pressure gradient
    [dpedr,dpedz] = gradient(pe(:,:,it).',dr,dz);
	dpedr = dpedr.';
    dpedz = dpedz.';

    % Plasma current density
    Jr(:,:,it+1)   = Jr(:,:,it) - dt*(Jr(:,:,it).*(v_ei + v_en + v_ion - v_loss) + e/me*(Bz.*Jphi(:,:,it) - Jz(:,:,it).*Bphi) - e/me*dpedr);
    Jphi(:,:,it+1) = Jphi(:,:,it) - dt*(Jphi(:,:,it).*(v_ei + v_en + v_ion - v_loss) + e/me*(Br.*Jz(:,:,it) - Jr(:,:,it).*Bz));
    Jz(:,:,it+1)   = Jz(:,:,it) - dt*(Jz(:,:,it).*(v_ei + v_en + v_ion - v_loss) + e/me*(Bphi.*Jr(:,:,it) - Jphi(:,:,it).*Br) - e/me*dpedz);
    
    %%
    %
    % 
    % Plasma density
    dn = n(:,:,it) + dt*(v_ion - v_loss);
    ind = find(dn<n0);
    dn(ind) = n0;
    n(:,:,it+1) = dn;
    n(1,:,it+1) = n0; n(end,:,it+1) = n0; n(:,1,it+1) = n0; n(:,end,it+1) = n0;

    % Electron pressure
    dn = pe(:,:,it) + dt*(2/3*(1 + (2*v_en + v_ion - v_loss)./v_ei/2).*rho.*(Jr(:,:,it).^2 + Jz(:,:,it).^2 + Jphi(:,:,it).^2) - 2*n(:,:,it)*e^2/mi.*rho.*(pe(:,:,it) - pion(:,:,it)));
    ind = find(dn<e*n0*Te0);
    dn(ind) = e*Te0*n0;
    pe(:,:,it+1) = dn;
    pe(1,:,it+1) = e*Te0*n0; pe(end,:,it+1) = e*Te0*n0; pe(:,1,it+1) = e*Te0*n0; pe(:,end,it+1) = e*Te0*n0;
    
    % Ion pressure
    dn = pion(:,:,it) + dt*(2*n(:,:,it)*e^2/mi.*rho.*(pe(:,:,it) - pion(:,:,it)));
    ind = find(dn<e*n0*Te0);
    dn(ind) = e*Te0*n0;
    pion(:,:,it+1) = dn;
        pion(1,:,it+1) = e*Te0*n0; pion(end,:,it+1) = e*Te0*n0; pion(:,1,it+1) = e*Te0*n0; pion(:,end,it+1) = e*Te0*n0;
    
    % Transport coefficients
    out = resistivity(n(:,:,it),pe(:,:,it)./n(:,:,it)/e,1);
    v_ei = out.v_ei;
    rho  = out.eta_par;
    v_ion = a.*sqrt(Jr(:,:,it).^2 + Jphi(:,:,it).^2 + Jz(:,:,it).^2)/e./n(:,:,it);
    v_loss = abs(Jr(:,:,it).*Br + Jphi(:,:,it).*Bphi + Jz(:,:,it).*Bz)./sqrt(Br.^2 + Bphi.^2 + Bz.^2)./n(:,:,it)/e./Leff;
%keyboard
end
keyboard
%% Plotting
figure(1)
clf
contourf(r,z,n(:,:,it+1).',51)
plot_tcabr(1)
hold off
colorbar

%%
if 0
clear all
syms Br Bz Bphi n e r Ephi
A = [e*n*r Bz -Bphi;-Bz e*n*r Br;Bphi -Br e*n*r];
J = simplify(inv(A)*[0; Ephi; 0]);
end

end