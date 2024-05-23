function breakdown_nova_pseudo

%% Input
clear all;
clc;
close all;
clf;
Ng=128;
Ntime =30*6;
dt = 1e-5; %tempo correto é 1e-4
R0 = 0.37;
B0 = 0.7; %aqui eu multipliquei por 3 no B0 para aumentar Bphi e por consequencia reduzir v_loss.
Vloop = 10;  %corrente gerando o campo E externo
p0 = 0.05; %0.01 % 1e-4 milibar = 0.01 Pa
ng = p0/1.38e-23/300; %particulas do gás
n0 = 1e8; %densidade de particulas inicial do plasma
Te0 = 0.026; % teperatura inicial do gas
gas = 'H2'; %tipo de gas
ligav_ion = 0;
% definindo a malha das coordenadas pseudotoridais
r0=0.07;
r_ini=0.001;
raio = linspace(r_ini,r0,Ng);
theta = linspace(0,2*pi,Ng); dtheta=2*pi/Ng; dR=(r0-r_ini)/Ng;
[X,A] = meshgrid(raio,theta); 
R = X.*cos(A)+R0;
Z = X.*sin(A); 
%
% Loading Green's functions table and plasma poloidal flux distribution
%load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65_pseudo.mat');
load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_128x128_pseudo.mat');
% Initializing physical constants and matrices
e = 1.6e-19;
me = 9.11e-31;
mi = 1.67e-27;

rg=R; zg=Z;
%dR = mean(diff(R(:,1))); 
dZ = mean(diff(Z(:,1))); %Definindo dr e dz
dr = dR; dz = dZ; %Definindo dr e dz

Br   = zeros(size(rg)); % definindo os Br e bz como zero
Bz = Br;

n   = n0*ones(size(rg,1),size(rg,2),Ntime+1); %definindo a condição inicial da densidade de particulas
eixo=[R0-0.07 R0+0.07 -0.07 0.07]; eixo2=[-0.07 +0.07 -0.07 0.07];
pe = e*n*Te0;  %a pressão de eletrons (pe) sai da densidade de particulas (n)  
pion = pe;     % a  pressão de ions (pion) começa igua a pressão de eletros
Jphi = n0*ones(size(rg,1),size(rg,2),Ntime+1); %definindo a densidade de corrente na direção toroidal como um valor inical constante bem pequeno 
Jz = Jphi.*0; % definindo a densidade de corrente na direção z como 0
Jr = Jphi.*0; % definindo a densidade de corrente na direção r como 0


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
% Cpl=I1/800;
% Iqpl = [Cpl*4 Cpl Cpl Cpl Cpl]*1e4; %correntes bobinas simulando plasma
% for icoil = 0:4
%     eval(['Br = Br + Ep' num2str(icoil) '.BR*Iqpl(' num2str(icoil+1) ');'])
%     eval(['Bz = Bz + Ep' num2str(icoil) '.BZ*Iqpl(' num2str(icoil+1) ');'])
% end

Bpol = sqrt(Br.^2 +Bz.^2);

if 1 figure(6); [~,h]=contourf(R,Z,Bpol,101); box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none');  ti1=['bpol(:,:)']; title(ti1);  end
%%
Bphi = R0*B0./rg;
Ephi = -Vloop/2/pi./rg;

% Computing quadrupole and toroidal magnetic and electric fields
out = resistivity_nova(n0,Te0,1);
a = first_townsend_coeff(p0,Ephi,gas,0);
Leff = R0*B0./(0.001+(Br.^2 + Bz.^2).^0.5);
rho = out.eta_par;
v_en = 7.89e11*p0;
v_ei = out.v_ei;
Jphi = Jphi.*(v_ei+v_en).*Ephi; 
%figure(4); [~,h]=contourf(R,Z,Jphi(:,:,1),101); hold on; axis equal; xlabel('R (m)'); ylabel('Z (m)'); colorbar;  axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); ti1=['jini(:,:,1)']; title(ti1); box on;
%keyboard
D = 0.3;  %termo de difusão
po=1;
count=1;

iuu=1;
s=10e8*exp(-100*((R-R0+0.02).^2+Z.^2)); %Aproximacao por uma Gaussiana
%figure(3);  colormap(jet); [q,h]=contourf(R,Z,s,101); hold on; axis equal; xlabel('R (m)'); ylabel('Z (m)'); colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); ti2=['s(:,:)']; title(ti2); 
%set(gcf, 'Position',  [1000, 1, 400, 400]); box on;
%keyboard
v_ion = ligav_ion*a.*sqrt(Jr(:,:,iuu).^2 + Jphi(:,:,iuu).^2 + Jz(:,:,iuu).^2)/e./n(:,:,iuu);
v_loss =ligav_ion*abs(Jr(:,:,iuu).*Br + Jphi(:,:,iuu).*Bphi + Jz(:,:,iuu).*Bz)./sqrt(Br.^2 + Bphi.^2 + Bz.^2)./n(:,:,iuu)/e./Leff;
%keyboard
for it = 1:Ntime
    %keyboard
    % Electron pressure gradient
    naux=size(pe(:,:,it));
    nR=naux(1);
    nZ=naux(2);
    dpedr=zeros(nR,nZ); dpedz=dpedr;
    for i=1:(nR-1)
    for j=1:(nZ-1)
       dpedr(i,j)=(pe(i+1,j,it)-pe(i,j,it))/dR/(raio(i));
       dpedz(i,j) = (pe(i,j+1,it)-pe(i,j,it))/(dtheta)/(raio(i));
    end
    end
    % Plasma current density
    interuptor1 = 0; %ativa/desatva   Jr*(v_ei... Se ativar interuptor 1 v ion e v loss explodem absurdamente rápido
    interuptor2 = 0; %ativa/desativa  e/me*(Bz.*Jphi(... ->se ativar int2 Jr e jz va oa infinito mt rapido ->uma vez que no calculo do intruptor 2 temos Bz-Br que no caso são 0, logo isso causa a explosao
    interuptor3 = 0; %ativa/desativa  e/me*dpedr ->se ativar int3 Jr e jz va oa infinito mt rapido
    %Jr(:,:,it+1)   = Jr(:,:,it) + dt*( - interuptor1*Jr(:,:,it).*(v_ei + v_en + v_ion - v_loss) + interuptor2*e/me*(Bz.*Jphi(:,:,it) - Jz(:,:,it).*Bphi) - interuptor3*e/me*dpedr);
    Jphi(:,:,it+1) = Jphi(:,:,it) + dt*(n(:,:,it)*e^2/me.*Ephi - interuptor1*Jphi(:,:,it).*(v_ei + v_en + v_ion - v_loss+s) + interuptor2*e/me*(Br.*Jz(:,:,it) - Jr(:,:,it).*Bz) );
    %Jphi(:,:,it+1) = -n(:,:,it)*e^2/me./(v_ei + v_en + v_ion - v_loss).*Ephi;
    %Jz(:,:,it+1)   = Jz(:,:,it) + dt*( - interuptor1*Jz(:,:,it).*(v_ei + v_en + v_ion - v_loss) + interuptor2*e/me*(Bphi.*Jr(:,:,it) - Jphi(:,:,it).*Br) - interuptor3*e/me*dpedz);
    
    %%
 
    % Plasma density
    laplaciano = D*Polar_laplaciano(n(:,:,it),dR,dtheta,raio);
    %keyboard
      laplacianoj = D*Polar_laplaciano(Jphi(:,:,it),dR,dtheta,raio);
    figure(8); [~,h]=contourf(X,A,laplacianoj);  box on;  axis equal; colorbar; axis(eixo2); set(h,'linecolor','none'); 
    dn = n(:,:,it) + dt*(v_ion - v_loss+laplaciano+s);
    %ind = find(dn<n0); dn(ind) = n0;
    n(:,:,it+1) = dn;
    n(1,:,it+1) = n(end,:,it+1); n(:,end,it+1) = n0;
        
    
    % Electron pressure
    auxx=2/3*(1 + (2*v_en + v_ion - v_loss+s)./v_ei/2).*rho.*(Jr(:,:,it).^2 + Jz(:,:,it).^2 + Jphi(:,:,it).^2);
    dn = pe(:,:,it) + dt*(0*auxx - 2*n(:,:,it)*e^2/mi.*rho.*(pe(:,:,it) - pion(:,:,it)));
    %ind = find(dn<e*n0*Te0); dn(ind) = e*Te0*n0;
    pe(:,:,it+1) = dn;
    pe(1,:,it+1) = pe(end,:,it+1); pe(:,end,it+1) = e*Te0*n0;
   
    % Ion pressure
    dn = pion(:,:,it) + dt*(2*n(:,:,it)*e^2/mi.*rho.*(pe(:,:,it) - pion(:,:,it)));
    %ind = find(dn<e*n0*Te0); dn(ind) = e*Te0*n0;
    pion(:,:,it+1) = dn;
    pion(1,:,it+1) = pion(end,:,it+1); 
    pion(:,end,it+1) = e*Te0*n0;
    % Transport coefficients
    iuu=it;
    %out = resistivity_nova(n(:,:,iuu),pe(:,:,iuu)./n(:,:,iuu)/e,1);
    %v_ei = out.v_ei;
  %  v_en = 7.89e11*(ng-n(:,:,iuu))*(1.38e-23)*300;   
   % rho  = out.eta_par;
    v_ion = ligav_ion*(a.*sqrt(Jr(:,:,iuu).^2 + Jphi(:,:,iuu).^2 + Jz(:,:,iuu).^2)/e./n(:,:,iuu));
    v_loss = ligav_ion*(abs(Jr(:,:,iuu).*Br + Jphi(:,:,iuu).*Bphi + Jz(:,:,iuu).*Bz)./sqrt(Br.^2 + Bphi.^2 + Bz.^2)./n(:,:,iuu)/e./Leff);
    %figure(5); [~,h]=contourf(R,Z,v_ion-v_loss,101); hold on; axis equal; xlabel('R (m)'); ylabel('Z (m)'); colorbar;  axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); ti1=['(v_ion-v_loss)(:,:,1)']; title(ti1);  
    %set(gcf, 'Position',  [1000, 401, 400, 400]); box on;
%     
%     max(max(v_ei)) 
%     max(max(-v_ei))
%     max(max(v_en)) 
%     max(max(-v_en))
  v_ion = ligav_ion*70000*v_ion; v_loss = v_loss*100*ligav_ion; 
  s=(10e8+10e11*(it/Ntime))*exp(-100*((R-0.37+0.02).^2+Z.^2));
    %keyboard
      if count == 3
          count=0;
    
          p1 = real(n(:,:,it));%./n(:,:,1);  %(0.00001+abs(n(:,:,1)));
          p2 = real(Jphi(:,:,it));%./Jphi(:,:,1); %(0.00001+abs(Jphi(:,:,1)));
          p3 = real(pe(:,:,it));%./pe(:,:,1); %(0.00001+abs(pe(:,:,1)));
          p4 = real(pion(:,:,it));%./pion(:,:,1); %(0.00001+abs(pion(:,:,1)));
         %disp(['v_ionz = ' num2str(max(max(v_ion))) ' particles/s']); disp(['v_loss = ' num2str(max(max(v_loss))) ' particles/s']); 
         %X Xu = [' it= ',num2str(it),' .']; disp(Xu); 
         figure(1);  set(gcf, 'Position',  [1, 1, 1200, 400]); box on; 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1,101); hold on; axis equal; xlabel('R (m)'); 
         ylabel('Z (m)'); colorbar; axis(eixo); set(h,'linecolor','none'); 
         ti1=['n(:,:,' num2str(it) ')0']; title(ti1); box on;
         subplot(1,2,2);  [~,h]=contourf(R,Z,p2,101); hold on; axis equal; xlabel('R (m)'); 
         ylabel('Z (m)'); colorbar;  axis(eixo); set(h,'linecolor','none'); 
         ti2=['Jphi(:,:,' num2str(it) ')0']; title(ti2);
         ti = [ti2 '.png']; saveas(gcf,ti)
         figure(2);  set(gcf, 'Position',  [1, 401, 1200, 400]); box on; 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p3,101); hold on; axis equal; xlabel('R (m)'); 
         ylabel('Z (m)'); colorbar;  axis(eixo); set(h,'linecolor','none'); 
         ti1=['pe(:,:,' num2str(it) ')0']; title(ti1); box on;
         subplot(1,2,2);  [~,h]=contourf(R,Z,p4,101); hold on; axis equal; xlabel('R (m)'); ylabel('Z (m)');
         colorbar;   axis(eixo); set(h,'linecolor','none'); 
         ti2=['pion(:,:,' num2str(it) ')0']; title(ti2);  
         ti = [ti2 '.png']; saveas(gcf,ti)
%          m1=max(max(-Jphi(:,:,it))) 
%          m2=max(max(pe(:,:,it))) 
         %m3=max(max(pion(:,:,it)))   %colormap(jet); [q,h]=contourf(R,Z,s+v_ion - v_loss,101); hold on; axis equal; xlabel('R (m)'); ylabel('Z (m)'); colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); ti2=['(v_i_o_n-v_l_o_s_s)(:,:,' num2str(it) ')']; title(ti2); 
        %F(it-1) = getframe;
         %pause(0.5);
       keyboard
       
     end
      count=count+1;
      disp(['s = ' num2str(max(max(s))) ' particles/s']);  disp(['v_ionz = ' num2str(max(max(v_ion))) ' particles/s']); disp(['v_loss = ' num2str(max(max(v_loss))) ' particles/s']); 
      Xu = [' it= ',num2str(it),' .']; disp(Xu); 
end

          p1 = n(:,:,:);
          p2 = real(Jphi(:,:,:));
          p3 = real(pe(:,:,:));
          p4 = real(pion(:,:,:));
          it = [Ntime/6 Ntime/3 3*Ntime/6 2*Ntime/3 5*Ntime/6 Ntime];
          
         figure(1);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(1)),101); 
         box on;  axis equal; colorbar; axis(eixoo); set(h,'linecolor','none'); 
         ti2=['ntod1']; ti1=['n(:,:,' num2str(it(1)) ')']; title(ti1); 
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(2)),101); 
         box on;  axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti = [ti2 '.png']; ti1=['n(:,:,' num2str(it(2)) ')']; title(ti1); saveas(gcf,ti)
         %%%%%%%%%%
         figure(2);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(3)),101); 
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti1=['n(:,:,' num2str(it(3)) ')']; title(ti1); 
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(4)),101);
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti2=['ntod2']; ti1=['n(:,:,' num2str(it(4)) ')']; title(ti1);
         ti = [ti2 '.png']; saveas(gcf,ti)
         %%%%%%%%%%
         figure(3);  set(gcf, 'Position',  [1, 401, 1200, 400]);  
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(5)),101);
         box on;  hold on; axis equal;  colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); %xlabel('R (m)'); ylabel('Z (m)');
         ti1=['n(:,:,' num2str(it(5)) ')']; title(ti1);
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(6)),101); 
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti1=['n(:,:,' num2str(it(6)) ')']; title(ti1);  ti2=['ntod3']; ti = [ti2 '.png'];  saveas(gcf,ti)
         
               p1=p2;    
         figure(1);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(1)),101); 
         box on;  axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti2=['jtod1']; ti1=['n(:,:,' num2str(it(1)) ')']; title(ti1); 
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(2)),101); 
         box on;  axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti = [ti2 '.png']; ti1=['Jphi(:,:,' num2str(it(2)) ')']; title(ti1); saveas(gcf,ti)
         %%%%%%%%%%
         figure(2);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(3)),101); 
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti1=['Jphi(:,:,' num2str(it(3)) ')']; title(ti1); 
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(4)),101);
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti2=['jtod2']; ti1=['n(:,:,' num2str(it(4)) ')']; title(ti1);
         ti = [ti2 '.png']; saveas(gcf,ti)
         %%%%%%%%%%
         figure(3);  set(gcf, 'Position',  [1, 401, 1200, 400]);  
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(5)),101);
         box on;  hold on; axis equal;  colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); %xlabel('R (m)'); ylabel('Z (m)');
         ti1=['Jphi(:,:,' num2str(it(5)) ')']; title(ti1);
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(6)),101); 
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti1=['Jphi(:,:,' num2str(it(6)) ')']; title(ti1);  ti2=['jtod3']; ti = [ti2 '.png'];  saveas(gcf,ti)
         
         p1=p3;
                   
         figure(1);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(1)),101); 
         box on;  axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti2=['petod1']; ti1=['pe(:,:,' num2str(it(1)) ')']; title(ti1); 
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(2)),101); 
         box on;  axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti = [ti2 '.png']; ti1=['pe(:,:,' num2str(it(2)) ')']; title(ti1); saveas(gcf,ti)
         %%%%%%%%%%
         figure(2);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(3)),101); 
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti1=['pe(:,:,' num2str(it(3)) ')']; title(ti1); 
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(4)),101);
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti2=['petod2']; ti1=['pe(:,:,' num2str(it(4)) ')']; title(ti1);
         ti = [ti2 '.png']; saveas(gcf,ti)
         %%%%%%%%%%
         figure(3);  set(gcf, 'Position',  [1, 401, 1200, 400]);  
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(5)),101);
         box on;  hold on; axis equal;  colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); %xlabel('R (m)'); ylabel('Z (m)');
         ti1=['n(:,:,' num2str(it(5)) ')']; title(ti1);
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(6)),101); 
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti1=['pe(:,:,' num2str(it(6)) ')']; title(ti1);  ti2=['petod3']; ti = [ti2 '.png'];  saveas(gcf,ti)
         
         p1=p4;
                   
         figure(1);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(1)),101); 
         box on;  axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti2=['pitod1']; ti1=['pion(:,:,' num2str(it(1)) ')']; title(ti1); 
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(2)),101); 
         box on;  axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti = [ti2 '.png']; ti1=['pion(:,:,' num2str(it(2)) ')']; title(ti1); saveas(gcf,ti)
         %%%%%%%%%%
         figure(2);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(3)),101); 
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti1=['pion(:,:,' num2str(it(3)) ')']; title(ti1); 
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(4)),101);
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti2=['pitod2']; ti1=['pion(:,:,' num2str(it(4)) ')']; title(ti1);
         ti = [ti2 '.png']; saveas(gcf,ti)
         %%%%%%%%%%
         figure(3);  set(gcf, 'Position',  [1, 401, 1200, 400]);  
         subplot(1,2,1);  [~,h]=contourf(R,Z,p1(:,:,it(5)),101);
         box on;  hold on; axis equal;  colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); %xlabel('R (m)'); ylabel('Z (m)');
         ti1=['pion(:,:,' num2str(it(5)) ')']; title(ti1);
         subplot(1,2,2);  [~,h]=contourf(R,Z,p1(:,:,it(6)),101); 
         box on;  hold on; axis equal; colorbar; axis([0.37-0.07 0.37+0.07 -0.07 0.07]); set(h,'linecolor','none'); 
         ti1=['pion(:,:,' num2str(it(6)) ')']; title(ti1);  ti2=['pitod3']; ti = [ti2 '.png'];  saveas(gcf,ti)
         
         
        % keyboard
end