close all                               % Fecha todos os gr�ficos
clear all                               % Exclui todas as vari�veis
clc                                     % Limpa a tela
format short                      % Formato para exibi��o num�rica
x = -1.05:0.1:1.05;
v = -1:0.1:1;
F = 0;
G = F;
[X,Y] = meshgrid(x);

for i = 1:21,
    p=i*0.1-1.1;
    F=F + 1./sqrt((X-p).^2+(Y+p).^2);
    G=G + 1./sqrt((X-p).^2+(Y-p).^2);
end
Fu = F+G+20./sqrt((X).^2+(Y).^2);;
[DX,DY] = gradient(Fu,.05,.05);
hold on
contour(X,Y,Fu)
quiver(X,Y,DX,DY)

X = -1:0.01:1;
y1 = sqrt(1-X.^2);
plot(X,y1,X,-y1)

for i= -1:0.1:1,
   scatter(i,sqrt(1-i^2))
   scatter(i,-sqrt(1-i^2))
end
