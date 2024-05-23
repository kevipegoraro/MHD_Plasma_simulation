function 
close all                               % Fecha todos os gr�ficos
clear all                               % Exclui todas as vari�veis
clc                                     % Limpa a tela
format short eng                        % Formato para exibi��o num�rica
x = -1.05:0.1:1.05;
y = -1.05:0.1:1.05;
[X,Y] = meshgrid(x,y);
F=1./(X+0.5);
F2=1./(X-0.5);
F4=1./(Y+0.5);
F5=1./(Y-0.5);
F3=20./sqrt((X).^2+(Y).^2);
[DX,DY] = gradient(F-F2+F3-F4-F5,.05,.05);
hold on
contour(X,Y,F-F2+F3+F4-F5)
quiver(X,Y,DX,DY)

X = -1:0.01:1;
y1 = sqrt(1-X.^2);
plot(X,y1,X,-y1)
plot(xx(1),yy(1),'ok')
plot(xx(2),-yy(2),'or')