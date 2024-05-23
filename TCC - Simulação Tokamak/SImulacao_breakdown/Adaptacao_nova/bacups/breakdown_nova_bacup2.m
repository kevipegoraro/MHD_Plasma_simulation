function breakdown_nova

%% Input
clear all
Ntime = 1000; %1000 ou 10000
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
load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65_2.mat');
%
% Initializing physical constants and matrices
e = 1.6e-19;
me = 9.11e-31;
mi = 1.67e-27;

[rg,zg] = meshgrid(r,z);  % criando a mash grid do r e z.
dr = mean(diff(r)); dz = mean(diff(z)); %Definindo dr e dz

Br   = zeros(size(rg)); % definindo os Br e bz como zero
Bz = Br;

n   = n0*ones(size(rg,1),size(rg,2),Ntime+1); %definindo a condição inicial da densidade de particulas

pe = e*n*Te0;  %a pressão de eletrons (pe) sai da densidade de particulas (n)  
pion = pe;     % a  pressão de ions (pion) começa igua a pressão de eletros

Jphi = 1e-6*ones(size(rg,1),size(rg,2),Ntime+1); %definindo a densidade de corrente na direção toroidal como um valor inical constante bem pequeno 
Jz = Jphi*0; % definindo a densidade de corrente na direção z como 0
Jr = Jphi*0; % definindo a densidade de corrente na direção r como 0


I1 = 1*100; %corrente passando pelas bobinas
fac = 0.4; % fator de mudança da corrente das bobinas verticias internas para as bobinas verticias externas.
% Current moments
Iq   = [-I1 I1 I1*fac -I1*fac]*1e4; % Quadrupole field


% calculando os campos das 4 bobinas verticais primarias 
for icoil = 1:4
    eval(['Br = Br + E' num2str(icoil) '.BR*Iq(' num2str(icoil) ');'])
    eval(['Bz = Bz + E' num2str(icoil) '.BZ*Iq(' num2str(icoil) ');'])
end

% calculando os campos das 4 bobinas que simulam o campo do plasma (isso é
% um aburdo fisico, pois anda não tenho o plasma mas vamo ve oq acontece)
% Cpl=-I1/800;
% Iqpl = [Cpl*4 Cpl Cpl Cpl Cpl]*1e4; %correntes bobinas simulando plasma
% for icoil = 0:4
%     eval(['Br = Br + Ep' num2str(icoil) '.BR*Iqpl(' num2str(icoil+1) ');'])
%     eval(['Bz = Bz + Ep' num2str(icoil) '.BZ*Iqpl(' num2str(icoil+1) ');'])
% end

Bpol = sqrt(Br.^2 +Bz.^2);

if 0
colormap(jet)
figure(2)
clf
[C,hh] = contourf(r,z,Bpol,101);
hold on
plot_nova(1)
hold off
set(hh,'linecolor','none')
axis([0.27 0.47 -0.1 0.1])
colorbar;
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
    interuptor1 = 0; %ativa/desatva   Jr*(v_ei...
    interuptor2 = 0; %ativa/desativa  e/me*(Bz.*Jphi(...
    interuptor3 = 0; %ativa/desativa  e/me*dpedr
    Jr(:,:,it+1)   = Jr(:,:,it) + dt*( - interuptor1*Jr(:,:,it).*(v_ei + v_en + v_ion - v_loss) + interuptor2*e/me*(Bz.*Jphi(:,:,it) - Jz(:,:,it).*Bphi) - interuptor3*e/me*dpedr);
    %Jphi(:,:,it+1) = Jphi(:,:,it) + dt*(n(:,:,it)*e^2/me.*Ephi - 0*Jphi(:,:,it).*(v_ei + v_en + v_ion - v_loss) + 0*e/me*(Br.*Jz(:,:,it) - Jr(:,:,it).*Bz) );
    Jphi(:,:,it+1) = -n(:,:,it)*e^2/me./(v_ei + v_en + v_ion - v_loss).*Ephi;
    Jz(:,:,it+1)   = Jz(:,:,it) + dt*( - interuptor1*Jz(:,:,it).*(v_ei + v_en + v_ion - v_loss) + interuptor2*e/me*(Bphi.*Jr(:,:,it) - Jphi(:,:,it).*Br) - interuptor3*e/me*dpedz);
    
    %%
 
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
    iuu=it;
    out = resistivity_nova(n(:,:,iuu),pe(:,:,iuu)./n(:,:,iuu)/e,1);
    v_ei = out.v_ei;
    v_en = 7.89e11*(ng-n(:,:,iuu))*(1.38e-23)*300;   
    rho  = out.eta_par;
    v_ion = 3*a.*sqrt(Jr(:,:,iuu).^2 + Jphi(:,:,iuu).^2 + Jz(:,:,iuu).^2)/e./n(:,:,iuu);
    v_loss = abs(Jr(:,:,iuu).*Br + Jphi(:,:,iuu).*Bphi + Jz(:,:,iuu).*Bz)./sqrt(Br.^2 + Bphi.^2 + Bz.^2)./n(:,:,iuu)/e./Leff;
end

%%
% for po = 2:round(Ntime/100,1):Ntime 
%     hold off
%     colormap(jet);figure(1);clf;
%     [q,h]=contourf(r,z,n(:,:,po),101);
%     hold on;
%     plot_nova(1);
%     colorbar
%     axis([0.27 0.47 -0.1 0.1]);
%     set(h,'linecolor','none'); 
%     pause(1.3); 
% end
disp(['v_ionz = ' num2str(max(max(v_ion))) ' particles/s'])
disp(['v_loss = ' num2str(max(max(v_loss))) ' particles/s'])
clf;
po=300;
subplot(2,2,1)
     colormap(jet);
     [q,h]=contourf(r,z,n(:,:,po),101);
     hold on;
     plot_nova(1);
     colorbar
     axis([0.27 0.47 -0.1 0.1]);
     set(h,'linecolor','none');
     title('n(:,:po)')

subplot(2,2,2)
     colormap(jet);
     [q,h]=contourf(r,z,Jphi(:,:,po),101);
     hold on;
     plot_nova(1);
     colorbar
     axis([0.27 0.47 -0.1 0.1]);
     set(h,'linecolor','none');
     title('Jphi(:,:po)')

subplot(2,2,3)
    colormap(jet);
     [q,h]=contourf(r,z,pion(:,:,po),101);
     hold on;
     plot_nova(1);
     colorbar
     axis([0.27 0.47 -0.1 0.1]);
     set(h,'linecolor','none');
     title('pion(:,:po)')

subplot(2,2,4)
    colormap(jet);
     [q,h]=contourf(r,z,pe(:,:,po),101);
     hold on;
     plot_nova(1);
     colorbar
     axis([0.27 0.47 -0.1 0.1]);
     set(h,'linecolor','none');
     titlee=num2str(po);
     legend(titlee);
     title('pe(:,:po)');
%text(0.15,0.45,num2str(po));
%%
keyboard
% colormap(jet);figure(1);clf;[q,h]=contourf(r,z,n(:,:,it+1),101);hold on;plot_nova(1);axis([0.27 0.47 -0.1 0.1]);set(h,'linecolor','none')
%figure(1);clf;[q,h]=contourf(r,z,Jphi(:,:,it+1),101);hold on;plot_nova(1);axis([0.27 0.47 -0.1 0.1]);set(h,'linecolor','none')
% %% Plotting

%%
if 0
clear all
syms Br Bz Bphi n e r Ephi
A = [e*n*r Bz -Bphi;-Bz e*n*r Br;Bphi -Br e*n*r];
J = simplify(inv(A)*[0; Ephi; 0]);
end

end