function plot_circulo3(corrente_plasma,corrente_bobinasparalelas1,corrente_bobinasparalelas2,corrente_bobinas)
close all                         
%clear all                              
%clc                                     
%format short                       
% exemplo de chamada: plot_circulo3(0.2,0.1/25,-0.1/20,0.05)
r=0.12/2; %(metros) raio da seção circular do tokamak nova furg 12 cm
h=r/100; %passo de plotagem
%escolhas das correntes
%corrente_plasma = 20; %potencia do campo da corrente de plasma
%corrente_bobinasparalelas = 10; %potencia do campo da corrente de plasma
%corrente_bobinas = 6;%potencia do campo da corrente de plasma
r1=2*r;
x = (-r1):h:(r1);
y = (-r1):h:(r1);
[X,Y] = meshgrid(x,y);

rd = r*(3/2-1/(6)); % rd deslocamento bobina direita
F2=0;
F1=0;
%4 campos devido as bobinas 
rb = r+r/4;
B1= corrente_bobinas./sqrt((X-rb).^2+(Y+rb).^2);
B2= corrente_bobinas./sqrt((X-rb).^2+(Y-rb).^2);
B3= corrente_bobinas./sqrt((X+rb).^2+(Y+rb).^2);
B4= corrente_bobinas./sqrt((X+rb).^2+(Y-rb).^2);
B_total =-B1+B2+B3-B4;
F_plasma=corrente_plasma./sqrt((X).^2+(Y).^2);
hold off

contour(X,Y,(F2-F1+F_plasma+B_total),400) %linhas de campo
hold on

%figure %abre nova figura
%plotagemde um circulo secção reta do Nova
alfa = 0:0.1:(2*pi);
xx = cos(alfa)*r;
yy = sin(alfa)*r;
plot(xx,yy,'r','linewidth',2)

plot([r*(1+1/3) r*(1+1/3)],[-2*r+0.1 2*r-0.1],'b','linewidth',2)
plot([-r*(3/2) -r*(3/2)],[-2*r+0.1 2*r-0.1],'b','linewidth',2)
%plot 4 bobinas
plot(-rb,rb,'xr','linewidth',3)
plot(-rb,-rb,'xr','linewidth',3)
plot(rb,rb,'xr','linewidth',3)
plot(rb,-rb,'xr','linewidth',3)
colorbar

axis([-(16/10)*r,(16/10)*r,-(16/10)*r,(16/10)*r]); %fixar o eixo
end