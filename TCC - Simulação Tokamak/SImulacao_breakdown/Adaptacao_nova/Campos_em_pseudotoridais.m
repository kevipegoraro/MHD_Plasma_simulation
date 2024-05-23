function Campos_em_pseudotoridais %(fig)
fig=1;


%close all %clear all
clc;
%criar variaveis
%%
load('Campos.mat');
%%
[R,Z] = meshgrid(r,z);
r0=0.12/2; p_central = 0.37; Ng=65; lines=Ng;
%---------
%criar vetores de posição em pseudotoridais
raio = linspace(0.001,2*r0,Ng);
theta = linspace(0,2*pi,Ng);
%cria a mashgrid psedoridal
[b11, b12] = meshgrid(raio,theta);
%------------------------

%plot do campos em cartesianas (R,Z considerando fhi = 0)
figure(fig)
colormap(jet); 
contour(r,z,Campos,lines); 
axis equal; 
plot_nova(fig); 
kevi = num2str(1);
text(0.43,0.45,kevi);
text(0,0.45,'ma hoi');
colorbar;
title('Campos magneticos poloidais em cartesianas')
%fazendo mudança de coordenadas de Campos
% lembrando q minha malha é r,z e fhi é 0.

%%%%% fazendo a figura em pseutoridais
%[b11, b12] = Pcoordenada(R,Z); %pega matrizes em coordenadas cartesianas e tranforma em matrizes em pseutoridais
[b21, b22] = Pcoordenada(Campos,Campos); 
figure(fig+1)
contour(b11, b12, Campos,lines); 
%quiver(b11, b12,b21, b22,0.6) %vetores
axis([0 2*r0 -pi pi]); %espaço de posições em pseutoridais
title('espaço de posições em pseutoridais')
figure(fig+2)
colormap(jet); 

contour(raio,theta,Campos,lines); 
%axis equal; 
plot_nova(fig); 
kevi = num2str(1);
text(0.43,0.45,kevi);
text(0,0.45,'ma hoi');
colorbar;

%%
%%% fazendo a figura voltando pra cartesianas
% [X1,Y1,Z1] = Pcoordenada_pseudo_to_cart(b11,b12,0); %pega matrezes de cordenadas pseutoridais e manda pra cartesianas
% [X2,Y2,Z2] = Pcoordenada_pseudo_to_cart(b21,b22,0); %pega matrezes de cordenadas pseutoridais e manda pra cartesianas
% figure(fig+2)
% %quiver(X1,Z1,X2,Z2,0.6);
% axis([p_central-2*r0 p_central+2*r0 -2*r0 2*r0]); %espaço de posições em cartesianas
% title('espaço de posições em pseutoridais trnformadas em cartesianas')
% 

end