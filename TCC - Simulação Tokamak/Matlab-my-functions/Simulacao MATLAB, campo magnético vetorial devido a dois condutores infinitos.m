%-------------------------------------------------------------------------%
%               UFSM - Universidade Federal de Santa Maria                %
%               Curso de Engenharia de Engenharia Elétrica                %
%                     ESP 1047 - Eletromagnetismo II                      %
%                                                                         %
%   Programadores:                                                        %
%       Laura Ferreira (201312207)                                        %
%                                                                         %
%                                                                         %
%   Versão: 1.0                                             30/09/2014    %
%=========================================================================%
%                        DESCRICAO DO PROGRAMA                            %
%                                                                         %
% Empregando o software MATLAB®, verificar o comportamento do campo       %
% magnetico vetorial devido à presença de dois condutores filamentares,   %
% infinitamente longos e paralelos.Calcular e exibir o vetor intensidadede%
% campo magnetico para um condutor retilineo.                             %
%=========================================================================%
%	Breve descrição do código.                                            %
%                                                                         %
%   v1.0 - Versão inicial.                                                %
%-------------------------------------------------------------------------%
close all                               % Fecha todos os gráficos
clear all                               % Exclui todas as variáveis
clc                                     % Limpa a tela
format short eng                        % Formato para exibição numérica
%-------------------------------------------------------------------------%
% ENTRADA DE DADOS
%-------------------------------------------------------------------------%
% Como entrada, o usuário deverá informar:
% 1. Os limites inferior e superior para o gráfico;
% 2. A distância de separação dos condutores;
% 3. o número de pontos para análise.
% 4. A intensidade da corrente elétrica em cada condutor e o respectivo sentido.
% Observar que as correntes podem assumir valores positivos (saindo do plano)
% ou negativo (entrando no plano);
% 
%-------------------------------------------------------------------------%

% _________________________________,,_____________________________________

%  1 LIMITES DE ANALISE

%  1.1 Eixo x

x_inf = input ('Insira o limite inferior para o grafico (m): ');
x_sup = input ('\nInsira o limite superior para o grafico (m): ');
x_lim = [x_inf x_sup]; %vetor limites

% Verificar se os dados foram inseridos

if isempty(x_lim)
    x_lim = [-10 10];
end

% 1.2 Eixo y (mesmos limites do eixo x)

y_lim = x_lim;

% _________________________________,,_____________________________________

% 2 DISTANCIA DE SEPARACAO DOS CONDUTORES

d=input('\nInsira a distancia de separacao dos condutores (m): ');

% Verificar se os dados foram inseridos

if isempty(d)
    d=5;
end

% _________________________________,,_____________________________________

% 3 NUMERO DE PONTOS PARA ANALISE

np=input('\n Insira o numero de pontos para analise: ');

% Verificar se os dados foram inseridos

if isempty(np)
    np=10;
end

% _________________________________,,_____________________________________

% 4 CORRENTE ELETRICA

fprintf('\n\nInsira a intensidade e o sentido das correntes eletricas nos condutores.');
fprintf('As correntes podem \nassumir valores positivos (saindo do plano) ou');
fprintf(' negativos entrando no plano.\n');


I1=input('\nCorrente no condutor 1: ');

% Verificar se os dados foram inseridos

if isempty(I1)
    I1=-4;
end

I2=input('\nCorrente no condutor 2: ');

% Verificar se os dados foram inseridos

if isempty(I2)
    I2=4;
end

%-------------------------------------------------------------------------%
% CÁLCULOS                                                                %
%-------------------------------------------------------------------------%

% Calculos em coordenadas cartesianas

% Cálculo dos pontos sobre os eixos x e y
eixo_x = linspace(x_lim(1), x_lim(2), np); %Matriz lim inf ao sup com np pontos
eixo_y = linspace(y_lim(1), y_lim(2), np);

% _________________________________,,_____________________________________

%                             VARREDURAS

for kk = 1:1:np                         % Varredura do eixo x
    for jj = 1:1:np                     % Varredura do eixo y
        
% _________________________________,,_____________________________________

        %CONDUTOR 1
        
        rho1 = [eixo_x(kk)+d/2 eixo_y(jj)];  % Vetor distância com relação ao eixo z
        
        modulo1 = norm(rho1);                % Modulo do vetor distancia
        
        a_phi1 = rho1/modulo1;               % Vetor unitario do distancia
        
        a_phi1 = [-a_phi1(2) a_phi1(1)];     % Rotacao a_phi = a_phi + 90º
        
  % _________________________________,,_____________________________________
        
        % CONDUTOR 2
          
        
        rho2 = [eixo_x(kk)-d/2 eixo_y(jj)];   % Vetor distância com relação ao eixo z
                                              % Está a "d" metros do condutor 1
                                          
        modulo2 = norm(rho2);                 % Modulo do vetor distancia
                         
        a_phi2 = rho2/modulo2;                % Vetor unitario do distancia
         
        a_phi2 = [-a_phi2(2) a_phi2(1)];      % Rotacao % a_phi = a_phi + 90º   
        
       
        H1 = I1/(2*pi*modulo1) * a_phi1;      % Vetor intensidade de campo magnético (A/m)
        H1x(jj,kk) = H1(1);                  % Componente x de H no ponto
        H1y(jj,kk) = H1(2);                  % Componente y de H no ponto
        
        H2 = I2/(2*pi*modulo2) * a_phi2;      % Vetor intensidade de campo magnético (A/m)
        H2x(jj,kk) = H2(1);                  % Componente x de H no ponto
        H2y(jj,kk) = H2(2);                  % Componente y de H no ponto
        
        Hx(jj,kk)=H1x(jj,kk)+H2x(jj,kk);   %Soma vetorial das componentes x de H1 e H2
        Hy(jj,kk)=H1y(jj,kk)+H2y(jj,kk);   %Soma vetorial das componentes y de H1 e H2
        
            end
end


% Cria matrizes com todas as combinações de posições
% X varia ao longo das colunas
% Y varia ao longo das linhas

[X, Y] = meshgrid(eixo_x, eixo_y);

%-------------------------------------------------------------------------%
% SAÍDA DE DADOS                                                          %
%-------------------------------------------------------------------------%

figure(1)
quiver(X, Y, Hx, Hy) %Plota a soma vetorial resultante
axis([x_lim y_lim])    %eixos iniciais

title('Comportamento do campo magnetico vetorial devido à presença de dois condutores filamentares');
% infinitamente longos e paralelos')
xlabel('Distância x (m)')
ylabel('Distância y (m)')


