function breakdown_nova

%% Input
clear all;
clc;
clf;
Ntime = 50; %1000 ou 10000
dt = 1e-4;
R0 = 0.37;
B0 = 0.7; %aqui eu multipliquei por 3 no B0 para aumentar Bphi e por consequencia reduzir v_loss.
Vloop = 10;  %corrente gerando o campo E externo
p0 = 0.05; %0.01 % 1e-4 milibar = 0.01 Pa
ng = p0/1.38e-23/300; %particulas do gás
n0 = 1e8; %densidade de particulas inicial do plasma
Te0 = 0.026; % teperatura inicial do gas
gas = 'H2'; %tipo de gas
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
Bz = Br; Bzpl = Br; Brpl = Br; 

n   = n0*ones(size(rg,1),size(rg,2),Ntime+1); %definindo a condição inicial da densidade de particulas

pe = e*n*Te0;  %a pressão de eletrons (pe) sai da densidade de particulas (n)  
pion = pe;     % a  pressão de ions (pion) começa igua a pressão de eletros

Jphi = 1e-6*ones(size(rg,1),size(rg,2),Ntime+1); %definindo a densidade de corrente na direção toroidal como um valor inical constante bem pequeno 
Jz = Jphi.*0; % definindo a densidade de corrente na direção z como 0
Jr = Jphi.*0; % definindo a densidade de corrente na direção r como 0


I1 = 1*100; %corrente passando pelas bobinas
fac = 0.345; % fator de mudança da corrente das bobinas verticias internas para as bobinas verticias externas.
% Current moments
Iq   = [-I1 I1 I1*fac -I1*fac]*1e4; % Quadrupole field


% calculando os campos das 4 bobinas verticais primarias 
for icoil = 1:4
    eval(['Br = Br + E' num2str(icoil) '.BR*Iq(' num2str(icoil) ');'])
    eval(['Bz = Bz + E' num2str(icoil) '.BZ*Iq(' num2str(icoil) ');'])
end

%calculando os campos das 4 bobinas que simulam o campo do plasma (isso é
%um aburdo fisico, pois anda não tenho o plasma mas vamo ve oq acontece)
Cpl=I1/800;
Iqpl = [Cpl*4 Cpl Cpl Cpl Cpl]*1e4; %correntes bobinas simulando plasma
for icoil = 0:4
    eval(['Brpl = Brpl + Ep' num2str(icoil) '.BR*Iqpl(' num2str(icoil+1) ');'])
    eval(['Bzpl = Bzpl + Ep' num2str(icoil) '.BZ*Iqpl(' num2str(icoil+1) ');'])
end

Bpol = sqrt(Br.^2 +Bz.^2); Bpol2 = sqrt(Brpl.^2 +Bzpl.^2); Bpol3 = sqrt(Br.^2 +Bz.^2+Brpl.^2 +Bzpl.^2);


colormap(jet)
figure(1)
clf
[C,hh] = contourf(r,z,Bpol,101);
hold on
plot_nova(1)
title('campo magnetico externo')
hold off
set(hh,'linecolor','none')
axis([0.27 0.47 -0.1 0.1])
colorbar;
colormap(jet)
figure(2)
clf
[C,hh] = contourf(r,z,Bpol2,101);
hold on
plot_nova(2)
title('campo magnetico plasma')
hold off
set(hh,'linecolor','none')
axis([0.27 0.47 -0.1 0.1])
colorbar;
colormap(jet)
figure(3)
clf
[C,hh] = contourf(r,z,Bpol3,101);
hold on
plot_nova(3)
title('campo magnetico externo mais campo plasma')
hold off
set(hh,'linecolor','none')
axis([0.27 0.47 -0.1 0.1])
colorbar;

%%

keyboard
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
D = 0.2;  %termo de difusão
po=1;
count=-1;
for it = 1:Ntime
    %keyboard
    % Electron pressure gradient
    [dpedr,dpedz] = gradient(pe(:,:,it),dr,dz);
	dpedr = dpedr;
    dpedz = dpedz;
    
    % Plasma current density
    interuptor1 = 0; %ativa/desatva   Jr*(v_ei... Se ativar interuptor 1 v ion e v loss explodem absurdamente rápido
    interuptor2 = 0; %ativa/desativa  e/me*(Bz.*Jphi(... ->se ativar int2 Jr e jz va oa infinito mt rapido ->uma vez que no calculo do intruptor 2 temos Bz-Br que no caso são 0, logo isso causa a explosao
    interuptor3 = 0; %ativa/desativa  e/me*dpedr ->se ativar int3 Jr e jz va oa infinito mt rapido
    %Jr(:,:,it+1)   = Jr(:,:,it) + dt*( - interuptor1*Jr(:,:,it).*(v_ei + v_en + v_ion - v_loss) + interuptor2*e/me*(Bz.*Jphi(:,:,it) - Jz(:,:,it).*Bphi) - interuptor3*e/me*dpedr);
    Jphi(:,:,it+1) = Jphi(:,:,it) + dt*(n(:,:,it)*e^2/me.*Ephi - interuptor1*Jphi(:,:,it).*(v_ei + v_en + v_ion - v_loss) + interuptor2*e/me*(Br.*Jz(:,:,it) - Jr(:,:,it).*Bz) );
    %Jphi(:,:,it+1) = -n(:,:,it)*e^2/me./(v_ei + v_en + v_ion - v_loss).*Ephi;
    %Jz(:,:,it+1)   = Jz(:,:,it) + dt*( - interuptor1*Jz(:,:,it).*(v_ei + v_en + v_ion - v_loss) + interuptor2*e/me*(Bphi.*Jr(:,:,it) - Jphi(:,:,it).*Br) - interuptor3*e/me*dpedz);
    
    %%
 
    % Plasma density
    dn = n(:,:,it) + dt*(v_ion - v_loss+D*del2(n(:,:,it),dr));
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
    v_ion = a.*sqrt(Jr(:,:,iuu).^2 + Jphi(:,:,iuu).^2 + Jz(:,:,iuu).^2)/e./n(:,:,iuu);
    v_loss = abs(Jr(:,:,iuu).*Br + Jphi(:,:,iuu).*Bphi + Jz(:,:,iuu).*Bz)./sqrt(Br.^2 + Bphi.^2 + Bz.^2)./n(:,:,iuu)/e./Leff;

    %keyboard
     if count == 1
         po = it; disp(['v_ionz = ' num2str(max(max(v_ion))) ' particles/s']); disp(['v_loss = ' num2str(max(max(v_loss))) ' particles/s']); count=0;
         Xu = [' it= ',num2str(po),' .']; disp(Xu)
         subplot(1,2,1); colormap(jet); [q,h]=contourf(r,z,n(:,:,po),101); hold on; plot_nova(1); colorbar ; axis([0.27 0.47 -0.1 0.1]); set(h,'linecolor','none'); title('n(:,:,po)');
         subplot(1,2,2); colormap(jet); [q,h]=contourf(r,z,Jphi(:,:,po),101); hold on; plot_nova(1); colorbar; axis([0.27 0.47 -0.1 0.1]); set(h,'linecolor','none'); title('Jphi(:,:,po)'); 
         %F(it-1) = getframe;
         keyboard
     end
      count=count+1;
end
%keyboard
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
po=2;
subplot(2,2,1)
     colormap(jet);
     [q,h]=contourf(r,z,n(:,:,po),101);
     hold on;
     plot_nova(1);
     colorbar
     axis([0.27 0.47 -0.1 0.1]);
     set(h,'linecolor','none');
     title('n(:,:,po)')

subplot(2,2,2)
     colormap(jet);
     [q,h]=contourf(r,z,Jphi(:,:,po),101);
     hold on;
     plot_nova(1);
     colorbar
     axis([0.27 0.47 -0.1 0.1]);
     set(h,'linecolor','none');
     title('Jphi(:,:,po)')

subplot(2,2,3)
    colormap(jet);
     [q,h]=contourf(r,z,pion(:,:,po),101);
     hold on;
     plot_nova(1);
     colorbar
     axis([0.27 0.47 -0.1 0.1]);
     set(h,'linecolor','none');
     title('pion(:,:,po)')

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
     title('pe(:,:,po)');
%text(0.15,0.45,num2str(po));
%%
%keyboard
% colormap(jet);figure(1);clf;[q,h]=contourf(r,z,n(:,:,it+1),101);hold on;plot_nova(1);axis([0.27 0.47 -0.1 0.1]);set(h,'linecolor','none')
%figure(1);clf;[q,h]=contourf(r,z,Jphi(:,:,it+1),101);hold on;plot_nova(1);axis([0.27 0.47 -0.1 0.1]);set(h,'linecolor','none')
% %% Plotting


end