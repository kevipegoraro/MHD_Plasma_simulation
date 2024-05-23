function plot_circulo(corrente_plasma,corrente_bobinasparalelas,corrente_bobinas,ativa_quiver)
%close all                               % Fecha todos os gr�ficos
%clear all                               % Exclui todas as vari�veis
%clc                                     % Limpa a tela
%format short                       % Formato para exibi��o num�rica
% exemplo plot_circulo(10,100,200,1) % plot_circulo(1,0.2,0.05,1)
r=0.12; %(metros) raio da seção circular do tokamak nova furg 12 cm
h=r/50; %passo de plotagem
%escolhas das correntes
%corrente_plasma = 20; %potencia do campo da corrente de plasma
%corrente_bobinasparalelas = 10; %potencia do campo da corrente de plasma
%corrente_bobinas = 6;%potencia do campo da corrente de plasma
x = (-2*r):h:(r*2);
y = (-2*r):h:(r*2);
[X,Y] = meshgrid(x,y);
F1=corrente_bobinasparalelas./(X+r*(3/2));
F2=corrente_bobinasparalelas./(X-r*(1+1/3));
%4 campos devido as bobinas 
rb = r+r/4;
B1= corrente_bobinas./sqrt((X-rb).^2+(Y+rb).^2);
B2= corrente_bobinas./sqrt((X-rb).^2+(Y-rb).^2);
B3= corrente_bobinas./sqrt((X+rb).^2+(Y+rb).^2);
B4= corrente_bobinas./sqrt((X+rb).^2+(Y-rb).^2);
B_total =-B1-B2-B3-B4;
F_plasma=corrente_plasma./sqrt((X).^2+(Y).^2);
hold off
if ativa_quiver == 1
    [DX,DY] = gradient((F2-F1+F_plasma+B_total),.24,.24);
    quiver(X,Y,DX,DY,15) %vetores
    hold on
end

contour(X,Y,(F2-F1+F_plasma+B_total),400) %linhas de campo
hold on

%figure %abre nova figura
%plotagemde um circulo secção reta do Nova
alfa = 0:0.1:(2*pi);
xx = cos(alfa)*r;
yy = sin(alfa)*r;
plot(xx,yy,'r','linewidth',2)
%plot bobinas paralelas
%reta = -2*r:0.1:4*r;
%plot(ones(length(reta)).*(r*(1+1/3)),reta,'b') %direita
%plot(ones(length(reta)).*(r*(1+1/3)+0.005),reta,'b')
plot([r*(1+1/3) r*(1+1/3)],[-2*r+0.1 2*r-0.1],'b','linewidth',2)
plot([-r*(3/2) -r*(3/2)],[-2*r+0.1 2*r-0.1],'b','linewidth',2)
%plot(ones(length(reta)).*(-r*(3/2)),reta,'b') %esqurda
%plot(ones(length(reta)).*(-r*(3/2)-0.005),reta,'b')

%plot 4 bobinas
plot(-rb,rb,'xr','linewidth',3)
plot(-rb,-rb,'xr','linewidth',3)
plot(rb,rb,'xr','linewidth',3)
plot(rb,-rb,'xr','linewidth',3)
colorbar
% 'or','ob','xr'
axis([-(16/10)*r,(16/10)*r,-(16/10)*r,(16/10)*r]); %fixar o eixo
end