function plot_circulo(corrente_plasma,corrente_bobinasparalelas1,corrente_bobinasparalelas2,corrente_bobinas)
close all                         
%clear all                              
%clc                                     
%format short                       
% exemplo de chamada: plot_circulo2(0.2,0.1/25,-0.1/20,0.05)
r=0.12/2; %(metros) raio da seção circular do tokamak nova furg 12 cm
h=r/100; %passo de plotagem
%escolhas das correntes
%corrente_plasma = 20; %potencia do campo da corrente de plasma
%corrente_bobinasparalelas = 10; %potencia do campo da corrente de plasma
%corrente_bobinas = 6;%potencia do campo da corrente de plasma

%F1=corrente_bobinasparalelas./(X+r*(3/2));
%F1 = [0 0 0 0 0 0 0 0 0 0 0];
%direita
r=0.12/2; 
h=r/100;
r1=2*r;
x = (-r1):h:(r1);
y = (-r1):h:(r1);
[X,Y] = meshgrid(x,y);
rd = r*(3/2-1/(6));
rb = r+r/4;% rd deslocamento bobina direita
% F_1=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(9/9)).^2);
% F_2=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(7/9)).^2);
% F_3=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(5/9)).^2);
% F_4=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(3/9)).^2);
% F_5=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(1/9)).^2);
%F_6=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(-1/9)).^2);
%F_7=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(-3/9)).^2);
%F_8=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(-5/9)).^2);
%F_9=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(-7/9)).^2);
%F_10=corrente_bobinasparalelas1./sqrt((X-rd).^2+(Y+r*(-9/9)).^2);
%F1=F_1+F_2+F_3+F_4+F_5+F_6+F_7+F_8+F_9+F_10;
%esquerda
%F2_1=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(9/9)).^2);
%F3_2=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(7/9)).^2);
%F4_3=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(5/9)).^2);
%F5_4=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(3/9)).^2);
%F6_5=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(1/9)).^2);
%F7_6=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(-1/9)).^2);
%%F8_7=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(-3/9)).^2);
%F9_8=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(-5/9)).^2);
%F10_9=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(-7/9)).^2);
%F11_10=corrente_bobinasparalelas2./sqrt((X+r*(3/2)).^2+(Y+r*(-9/9)).^2);
%F2=F2_1+F3_2+F4_3+F5_4+F6_5+F7_6+F8_7+F9_8+F10_9+F11_10;
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

plot(-rb,rb,'xr','linewidth',3)
plot(-rb,-rb,'xr','linewidth',3)
plot(rb,rb,'xr','linewidth',3)
plot(rb,-rb,'xr','linewidth',3)
colorbar
axis([-(16/10)*r,(16/10)*r,-(16/10)*r,(16/10)*r]); %fixar o eixo
end