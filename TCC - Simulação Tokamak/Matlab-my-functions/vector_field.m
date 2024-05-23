
close all                               % Fecha todos os gr�ficos
clear all                               % Exclui todas as vari�veis
clc                                     % Limpa a tela
format short 
l=-1.05:0.1:1.05; %vetor principal

[X,Y] = meshgrid(l); %grid para plotar circulo
Z = 1./sqrt(X.^2 + Y.^2);


[x,y] = meshgrid(l);  %grid para plotar o campo vetorial
u=-cos(x).*sin(y)./sqrt(X.^2 + Y.^2)
v=cos(y).*sin(x)./sqrt(X.^2 + Y.^2)
hold on

quiver(x,y,u,v)
contour(x,y,Z)

%plotar um circulo:
X = -1:0.01:1;
y1 = sqrt(1-X.^2);
plot(X,y1,X,-y1)