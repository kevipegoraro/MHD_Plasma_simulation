%Entendendo Vetores em coordenadas Pseudoridais

%oque eu sei: cada vetor E tera um formato E=E_r+E_theta+E_fhi
%todo vetor começa em um ponto A_i(r1,theta1,fhi1) e termina em A_f(r2,theta2,fhi2) 
%onde a norma(E) é a intensidade do campo no ponto.
%O vetor E que começa em A_i e termina em A_f é dado por: E =
%vv{r}*(r2-r1)+vv{\theta}*(theta2-theta1)+vv{\fhi}*(fhi2-fhi1)
%Cada campo elétrico começa no ponto que estou olhando seu valor, e termina
%no ponto que determina sua direção e módulos que ja são conhecidos.
%assim posso desenhar via A_i e A_f os vetores do meu campo elétrico
%para calcular a dintancia entre 2 pontos em cordenadas pseudotoridais irei
%passar os dois pontos para cartesianas e então fazer o calculo da norma da subtraçã ode ambos os vetores.
% v=[0.03 0 0]; %vetor começa em v e termina em l 
% l=[0.06 pi pi/2];
% X1=pseudotoridal2(v)
% X2=pseudotoridal2(l)
% Norma=norm(X2-X1) 

clear all;
clc;
clf;
Ntime = 50; %1000 ou 10000
dt = 1e-4;
R0 = 0.37;
B0 = 3*0.7; %aqui eu multipliquei por 3 no B0 para aumentar Bphi e por consequencia reduzir v_loss.
Vloop = 10;  %corrente gerando o campo E externo
p0 = 0.05; %0.01 % 1e-4 milibar = 0.01 Pa
ng = p0/1.38e-23/300; %particulas do gás
n0 = 1e8; %densidade de particulas inicial do plasma
Te0 = 0.026; % teperatura inicial do gas
gas = 'H2'; %tipo de gas

Ng=65; r0=0.07; 
fhi = linspace(0,2*pi,Ng);
Mo = pi*4E-7;
corrente = Mo*4/0.001;
p_central = 0.37;
%fazendo ponto a ponto as equipontenciais de campo em pseudotoridais
r_ini=0.001;
raio = linspace(r_ini,r0,Ng);
theta = linspace(0,2*pi,Ng);
[R,A] = meshgrid(raio,theta);
%E=corrente./sqrt(((0.37+R*cos(A))).^2+(R*sin(A)).^2);
x = R.*cos(A)+0.37;
y = R.*sin(A); 
%
% Loading Green's functions table and plasma poloidal flux distribution
load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65_pseudo.mat');
%
% Initializing physical constants and matrices
e = 1.6e-19;
me = 9.11e-31;
mi = 1.67e-27;


dr = mean(diff(x)); dz = mean(diff(y)); %Definindo dr e dz

Br   = zeros(size(x)); % definindo os Br e bz como zero
Bz = Br;

I1 = 1*100; %corrente passando pelas bobinas
fac = 0.4; % fator de mudança da corrente das bobinas verticias internas para as bobinas verticias externas.
% Current moments
Iq   = [-I1 I1 I1*fac -I1*fac]*1e4; % Quadrupole field


% calculando os campos das 4 bobinas verticais primarias 
for icoil = 1:4
    eval(['Br = Br + E' num2str(icoil) '.BR*Iq(' num2str(icoil) ');'])
    eval(['Bz = Bz + E' num2str(icoil) '.BZ*Iq(' num2str(icoil) ');'])
end
Bpol = sqrt(Br.^2 +Bz.^2);

if 1
colormap(jet); [q,h]=contourf(x,y,Bpol,101); 
colorbar ;  
set(h,'linecolor','none'); 
title('campo magnetico poloidalr');
end
%%
Bphi = R0*B0./x;
Ephi = -Vloop/2/pi./x;

% 
% figure(2)
% E=corrente./sqrt((x+0.37).^2+(y).^2);
% colormap(jet); [q,h]=contourf(x,y,E,101); 
% hold on; 
% colorbar ;  
% set(h,'linecolor','none'); 
% title('plot contuo em polar');
% 
% figure(3)
% colormap(jet); [q,h]=contourf(x,y,Ephi,101); 
% hold on; 
% colorbar ;  
% set(h,'linecolor','none'); 
% title('plot contuor campo eletrico externo');
% 
% 
% figure(4)
% colormap(jet); [q,h]=contourf(x,y,Bphi,101); 
% hold on; 
% colorbar ;  
% set(h,'linecolor','none'); 
% title('plot contuor campo magnetico externo');







    %fronteira no eixo x e y
%     for lo=1:1:Ng
%         for li=1:1:Ng
%             if  raio(lo)*cos(theta(li))+0.37==zg(lo,li)
%                 n(lo,li,it+1)=n0;
%             end
%             if  raio(lo)*sin(theta(li))==rg(lo,li)
%                 n(lo,li,it+1)=n0;
%             end
%         end
%     end






















% x = linspace(p_central-r0,p_central+r0,Ng);
% % y = linspace(-r0,r0,Ng);
% 
% for i=1:1:Ng
%     for j=1:1:Ng
%          v=[raio(i) theta(j) 0]; %ponto q estou calculando o campo
%          l=[R_d pi 0]; %localização em pseudotoridais do campo pontual
%          %Ii=pseudotoridal2(v); Ff=pseudotoridal2(l); x_i=Ii(1); y_i=Ii(2);
%          %x_f=R_d; %Ff(1); y_f=0; quad = atan((y_f-y_i)/(x_f-x_i));
%          %Raios(i,j) = abs(x_f-x_i)/cos(quad); %
%          E(i,j)=corrente*(1./tocartotonorma(v,l)^2); %calcuo do vetor campo em si
%          %Angulos(i,j) =
%          %pi-quad;%atan(raio(i)*sin(pi-theta(j))/(R_d-raio(i)*cos(pi-theta(j))));
%          % mesmo campo e em cartesianas
% %          v2=[x(i) y(j) 0]; %ponto q estou calculando o campo
% %          l2=[p_central+0.02 0 0]; %localização em cartesianas do campo
% %          pontual E2(i,j) = corrente*(1./norm(v2-l2)^2); %calcuo do vetor
% %          campo em si
%    end
% end

%z=R-A;

% Hide the POLAR function data and leave annotations

% Turn off axes and set square aspect ratio
% [DX, DY] = gradient(Angulos, 0.01, 0.01);
% hold on
% %contour(raio, theta, E,Ng)
% plot(0.04*ones(Ng),theta)
% plot((r_ini+0.01)*ones(Ng),theta)
% plot(raio,pi/2*ones(Ng))
% plot(raio,pi*ones(Ng))
% colormap(jet);
% colorbar;
% % caxis([1 9]*10^3)
% %quiver(R, A, Raios,Angulos,0.01)
% quiver(R, A, DX, DY)
% % %hold on;  grid off;
% % figure(2);
% % 
% % contour(x, y, E,Ng)

%%

% %fazendo ponto a ponto as equipontenciais de campo em cartecianas
% x = linspace(p_central-t_contour*r0,p_central+t_contour*r0,Ng);
% y = linspace(p_central-t_contour*r0,p_central+t_contour*r0,Ng);
% for i=1:1:Ng
%     for j=1:1:Ng
%          E(i,j) = corrente*(1./sqrt((x(i)-p_central).^2+y(j).^2));
%     end
% end
% %%
% [DX, DY] = gradient(E, 0.1, 0.1);
% contour(x, y, E)
% hold on;  grid off;
% quiver(x, y, DX, DY)
% %%


