function breakdown_nova

%% Input
clear all
Ntime = 10000;
dt = 1e-5;
R0 = 0.37;
B0 = 0.7;
Vloop = 10;
p0 = 0.05; %0.01 % 1e-4 milibar = 0.01 Pa
ng = p0/1.38e-23/300;
n0 = 1e8;
Te0 = 0.026;
gas = 'H2';
%
% Loading Green's functions table and plasma poloidal flux distribution
%load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65.mat');
load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65_2.mat');
%
% Initializing physical constants and matrices
e = 1.6e-19;
me = 9.11e-31;
mi = 1.67e-27;
[rg,zg] = meshgrid(r,z); dr = mean(diff(r)); dz = mean(diff(z));
Br   = zeros(size(rg)); Bz = Br;
n   = n0*ones(size(rg,1),size(rg,2),Ntime+1); pe = e*n*Te0; pion = pe;
Jphi = 1e-6*ones(size(rg,1),size(rg,2),Ntime+1); Jz = Jphi*0; Jr = Jphi*0;


I1 = 1*100;
fac = 0.4;
% Current moments
Iq   = [-I1 I1 I1*fac -I1*fac]*1e4; % Quadrupole field
%Iq   = [1 1 1 1 1 1 -1 -1 -0.6 -0.6 0.6 0.6 1]*1e4;  % Quadrupole field

% Computing quadrupole and toroidal magnetic and electric fields
for icoil = 1:4
    eval(['Br = Br + E' num2str(icoil) '.BR*Iq(' num2str(icoil) ');'])
    eval(['Bz = Bz + E' num2str(icoil) '.BZ*Iq(' num2str(icoil) ');'])
end
% for icoil = 5:10
%     eval(['Br = Br + E' num2str(icoil) '.BR*Iq(' num2str(icoil) ');'])
%     eval(['Bz = Bz + E' num2str(icoil) '.BZ*Iq(' num2str(icoil) ');'])
% end
Bpol = sqrt(Br.^2 +Bz.^2);

if 0
colormap(jet)
figure(1)
clf
[C,hh] = contourf(r,z,Bpol,101);
hold on
plot_nova(1)
hold off
set(hh,'linecolor','none')
axis([0.27 0.47 -0.1 0.1])
end
%%
Bphi = R0*B0./rg;
Ephi = -Vloop/2/pi./rg;

% Computing quadrupole and toroidal magnetic and electric fields
out = resistivity_nova(n0,Te0,1);
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
    [dpedr,dpedz] = gradient(pe(:,:,it),dr,dz);
	dpedr = dpedr;
    dpedz = dpedz;
    
    % Plasma current density
    Jr(:,:,it+1)   = Jr(:,:,it) + dt*( - 0*Jr(:,:,it).*(v_ei + v_en + v_ion - v_loss) + 0*e/me*(Bz.*Jphi(:,:,it) - Jz(:,:,it).*Bphi) - 0*e/me*dpedr);
    %Jphi(:,:,it+1) = Jphi(:,:,it) + dt*(n(:,:,it)*e^2/me.*Ephi - 0*Jphi(:,:,it).*(v_ei + v_en + v_ion - v_loss) + 0*e/me*(Br.*Jz(:,:,it) - Jr(:,:,it).*Bz));
    Jphi(:,:,it+1) = -n(:,:,it)*e^2/me./(v_ei + v_en + v_ion - v_loss).*Ephi;
    Jz(:,:,it+1)   = Jz(:,:,it) + dt*( - 0*Jz(:,:,it).*(v_ei + v_en + v_ion - v_loss) + 0*e/me*(Bphi.*Jr(:,:,it) - Jphi(:,:,it).*Br) - 0*e/me*dpedz);
    
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
    out = resistivity_nova(n(:,:,it+1),pe(:,:,it+1)./n(:,:,it+1)/e,1);
    v_ei = out.v_ei;
    v_en = 7.89e11*(ng-n(:,:,it+1))*(1.38e-23)*300;
    rho  = out.eta_par;
    v_ion = a.*sqrt(Jr(:,:,it+1).^2 + Jphi(:,:,it+1).^2 + Jz(:,:,it+1).^2)/e./n(:,:,it+1);
    v_loss = abs(Jr(:,:,it+1).*Br + Jphi(:,:,it+1).*Bphi + Jz(:,:,it+1).*Bz)./sqrt(Br.^2 + Bphi.^2 + Bz.^2)./n(:,:,it+1)/e./Leff;
end
disp(['v_ionz = ' num2str(max(max(v_ion))) ' particles/s'])
disp(['v_loss = ' num2str(max(max(v_loss))) ' particles/s'])
keyboard

% %% Plotting
% figure(1)
% clf
% i=4;
% contourf(r,z,n(:,:,i).',51)
% color bar;
% 
% figure(2)
% %contourf(r,z,Jr(:,:,i).',51)
% contourf(r,z,Jr(:,:,i),51)
% figure(3)
% %contourf(r,z,Jphi(:,:,i).',51)
% contourf(r,z,Jphi(:,:,i),51)
% figure(4)
% contourf(r,z,Jz(:,:,i).',51)
% figure(5)
% contourf(r,z,sqrt(Jz(:,:,i).^2+Jz(:,:,i).^2+Jphi(:,:,i).^2).',51)
% 
% %contourf(r,z,n(:,:,it+1).',51)
% plot_nova(1)
% plot_nova(2)
% plot_nova(3)
% plot_nova(4)
% plot_nova(5)
% %colorbar

%%
if 0
clear all
syms Br Bz Bphi n e r Ephi
A = [e*n*r Bz -Bphi;-Bz e*n*r Br;Bphi -Br e*n*r];
J = simplify(inv(A)*[0; Ephi; 0]);
end

end