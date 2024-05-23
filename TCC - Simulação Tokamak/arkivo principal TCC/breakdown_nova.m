function breakdown_nova
% Codigo explicito no tempo
% Autoria Professor Gustavo P Canal - Adaptado para o NOVA neste trabalho
Ntime = 50; dt = 1e-6; % Numero de tempos para simular e intervalo de tempo
R0 = 0.37; B0 = 0.7; p0 = 0.05; % Dados NOVA
Vloop = 10; % Corrente gerando o campo E externo
ng = p0/1.38e-23/300; % Particulas do gas
n0 = 1e8; % Densidade de particulas inicial do plasma
Te0 = 0.026; % Temperatura inicial do gas
gas = 'H2'; % Tipo de gas
%
% Carregando as funcoes de Grren e fluxo de plasma poloidal
load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65_2.mat');
%
% Definindo as constantes fisicas e matrizes
e = 1.6e-19;
me = 9.11e-31;
mi = 1.67e-27;

[rg,zg] = meshgrid(r,z);  % Criando a meshgrid do r e z.
dr = mean(diff(r)); dz = mean(diff(z)); % Definindo dr e dz

Br   = zeros(size(rg)); % Definindo os Br e Bz como zero
Bz = Br;

n   = n0*ones(size(rg,1),size(rg,2),Ntime+1); % Definindo a condição inicial da densidade de particulas

pe = e*n*Te0;  % A pressao de eletrons (pe) sai da densidade de particulas (n)
pion = pe;     % A pressão de ions (pion) começa igua a pressão de eletros

Jphi = 1e-6*ones(size(rg,1),size(rg,2),Ntime+1); % Definindo a densidade de corrente na direção toroidal como um valor inical constante bem pequeno
Jz = Jphi.*0; % Definindo a densidade de corrente na direcao z como 0
Jr = Jphi.*0; % Definindo a densidade de corrente na direcao r como 0
corrente=0; it=1; Bpl_r=ones(size(rg,1),size(rg,2),Ntime+1); Bpl_z=Bpl_r; Bpl_phi=Bpl_r;
for ii=1:10
    corrente=corrente+Jphi(26+ii,26+ii,1);
end
corrente=corrente*(dr*10)^3*e*2;
Bpl_r(:,:,it) = Ep0.BR*corrente;
Bpl_z(:,:,it) = Ep0.BZ*corrente;
Bpl_phi(:,:,it)=Ep0.G*corrente*0;

I1 = 1*100; % Corrente passando pelas bobinas
fac = 0.4; % Fator de mudança da corrente das bobinas verticias internas para as bobinas verticias externas.
% Momentos de corrente
Iq   = [-I1 I1 I1*fac -I1*fac]*1e4;


% Calculando os campos das 4 bobinas verticais primarias
for icoil = 1:4
    eval(['Br = Br + E' num2str(icoil) '.BR*Iq(' num2str(icoil) ');'])
    eval(['Bz = Bz + E' num2str(icoil) '.BZ*Iq(' num2str(icoil) ');'])
end

Bphi = R0*B0./rg;
Ephi = -Vloop/2/pi./rg;
% Iniciando os coeficientes de transporte
out = resistivity_nova(n0,Te0,1);
a = first_townsend_coeff(p0,Ephi,gas,0);
Leff = R0*B0./(0.001+(Br.^2 + Bz.^2).^0.5);
rho = out.eta_par;
v_en = 7.89e5*p0;
v_ion = 0.01;
v_loss = 0.005;
v_ei = out.v_ei;
D = 0.2;  %termo de difusao
po=1;
count=-2;
for it = 1:Ntime
    
    % Gradiente de pressao eletronica
    [dpedr,dpedz] = gradient(pe(:,:,it),dr,dz);
    dpedr = dpedr;
    dpedz = dpedz;
    
    % Densidade de corrente de plasma
    modJ=sqrt(Jr(:,:,it).^2 + Jphi(:,:,it).^2 + Jz(:,:,it).^2);
    Jr(:,:,it+1)   = 0;
    Jz(:,:,it+1)   = 0;
    part1=-Jphi(:,:,it).*(v_ei  + v_ion - v_loss);
    part2=e/me*((Br+Bpl_r(:,:,it)).*Jz(:,:,it) - Jr(:,:,it).*(Bz+Bpl_z(:,:,it)));
    part3=n(:,:,it)*e^2/me.*Ephi;
    Jphi(:,:,it+1) = Jphi(:,:,it) + dt*(part1+part2+part3);
    % Uma aproximacao dos campos magnetico e eletrico gerados pela corrente
    % de plasma em cada tempo
    corrente=0; % Atualizando os campos magneticos do plasma
    for ii=1:10
        corrente=corrente+n(26+ii,26+ii,it);
    end
    corrente=corrente*e*2;
    Bpl_r(:,:,it) = Ep0.BR*corrente;
    Bpl_z(:,:,it) = Ep0.BZ*corrente;
    Bpl_phi(:,:,it)=Ep0.G*corrente;
    
    
    % Densidade de particulas
    div=calcdiv(Jr(:,:,it),Jz(:,:,it),dr,dz);
    dn = n(:,:,it) + dt*(div+v_ion - v_loss+D*del2(n(:,:,it),dr));
    ind = find(dn<n0);
    dn(ind) = n0;
    n(:,:,it+1) = dn;
    n(1,:,it+1) = n0; n(end,:,it+1) = n0; n(:,1,it+1) = n0; n(:,end,it+1) = n0;
    
    % Pressao eletronica
    dn = pe(:,:,it) + dt*(2/3*(1 + (2*v_en + v_ion - v_loss)./v_ei/2).*rho.*(modJ.^2) - 2*n(:,:,it)*e^2/mi.*rho.*(pe(:,:,it) - pion(:,:,it)));
    ind = find(dn<e*n0*Te0);
    dn(ind) = e*Te0*n0;
    pe(:,:,it+1) = dn;
    pe(1,:,it+1) = e*Te0*n0; pe(end,:,it+1) = e*Te0*n0; pe(:,1,it+1) = e*Te0*n0; pe(:,end,it+1) = e*Te0*n0;
    
    % Pressao ionica
    dn = pion(:,:,it) + dt*(2*n(:,:,it)*e^2/mi.*rho.*(pe(:,:,it) - pion(:,:,it)));
    ind = find(dn<e*n0*Te0);
    dn(ind) = e*Te0*n0;
    pion(:,:,it+1) = dn;
    pion(1,:,it+1) = e*Te0*n0; pion(end,:,it+1) = e*Te0*n0; pion(:,1,it+1) = e*Te0*n0; pion(:,end,it+1) = e*Te0*n0;
    
    % Coeficientes de transporte
    out = resistivity_nova(n(:,:,it),pe(:,:,it)./n(:,:,it)/e,1);
    v_ei = out.v_ei;
    v_en = 7.89e11*(ng-n(:,:,it))*(1.38e-23)*300;
    rho  = out.eta_par;
    v_ion = a.*modJ/e./n(:,:,it);
    v_loss = abs(Jr(:,:,it).*(Br+Bpl_r(:,:,it)) + Jphi(:,:,it).*(Bphi+Bpl_phi(:,:,it)) + Jz(:,:,it).*Bz)./sqrt((Br+Bpl_r(:,:,it)).^2 + (Bphi+Bpl_phi(:,:,it)).^2 + (Bz+Bpl_z(:,:,it)).^2)./n(:,:,it)/e./Leff;
    
    % Plotagem
    if count == 1
        po = it; disp(['v_ionz = ' num2str(max(max(v_ion))) ' particles/s']); disp(['v_loss = ' num2str(max(max(v_loss))) ' particles/s']); count=0;
        Xu = [' it= ',num2str(po),' .']; disp(Xu)
        subplot(2,1,1); colormap(jet); [q,h]=contourf(r,z,n(:,:,po),101); hold on; plot_nova(1); colorbar ; axis([0.27 0.47 -0.1 0.1]); set(h,'linecolor','none');  Xu1=['n(:,:,' num2str(po) ')']; title(Xu1);
        subplot(2,1,2); colormap(jet); [q,h]=contourf(r,z,Jphi(:,:,po),101); hold on; plot_nova(1); colorbar; axis([0.27 0.47 -0.1 0.1]); set(h,'linecolor','none'); Xu2=['Jphi(:,:,' num2str(po) ')']; title(Xu2);
        hold off
        figure(4);  colormap(jet); [q,h]=contourf(r,z,sqrt((Br+Bpl_r(:,:,it)).^2 + (Bz+Bpl_z(:,:,it)).^2),101); colorbar ; axis([0.27 0.47 -0.1 0.1]); set(h,'linecolor','none');  Xu1=['B_pol_total']; title(Xu1); axis equal;
        figure(5);  colormap(jet); [q,h]=contourf(r,z,sqrt((Bphi+Bpl_phi(:,:,it)).^2),101); colorbar; axis([0.27 0.47 -0.1 0.1]); set(h,'linecolor','none');  Xu1=['B_toroidal_total']; title(Xu1); axis equal;
        keyboard
    end
    count=count+1;
end
end