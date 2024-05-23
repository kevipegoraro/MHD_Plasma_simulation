close all                               % Fecha todos os gr�ficos
clear all                               % Exclui todas as vari�veis
clc                                     % Limpa a tela
format short eng                        % Formato para exibi��o num�rica
x = -1.05:0.1:1.05;
z = x;
y = -1.05:0.1:2.05;
[X,Y,Z] = meshgrid(x,y,z);
F=10./sqrt(X.^2+Y.^2+Z.^2);
F2=10./sqrt(X.^2+(Y-1).^2+Z.^2);
F3=10./sqrt((X+0.5).^2+(Y+0.2).^2+(Z-0.2).^2);
[DX,DY,DZ] = gradient(F+F2-F3,.05,.05,.05);
quiver3(X,Y,Z,DX,DY,DZ)

%F = X.*exp(-X.^2-Y.^2);
%surf(X,Y,F)