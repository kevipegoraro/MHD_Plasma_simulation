function Simulacao_n_J_Gaussiana
clear all
%close all
clc
global G e me mi CAMPO
Neq=4; e= 1.6e-19; me = 9.11e-31; mi = 1.67e-27;
Ntime = 50; 
dt = 1e-7; tlist = 0:dt:dt*Ntime;
G.R0 = 0.37;  G.control=0;  G.dt=dt;
G.B0 = 0.7;   G.p0 = 0.05; 
G.ng = G.p0/1.38e-23/300;  G.n0 = 1e8;
G.Te0 = 0.026;             G.gas = 'H2'; 
G.Vloop = 10*3;              G.a0=0.06;  G.eixo=[G.R0-G.a0 G.R0+G.a0 -G.a0 G.a0];
G.C1 = 1;                   G.V = 1;     G.camp=1;  
G.ligacorrente = 1; % Ativa o calculo do campo magnetico gerado pelo plasma
G.cplott = 0; G.plotton = 1; % Ativa o plot da fonte em todos os tempos
%%
% u(1,:) -> n
% u(2,:) -> J_phi
% u(3,:) -> pe
% u(4,:) -> pion
model = createpde(Neq);  % dl = decsg(gd,sf,ns)  pdecirc(0.37,0,0.06)  pdecirc(0.37,0,0.01)
dl=[1 1 1 1; 0.31 0.37 0.43 0.37; 0.37 0.43 0.37 0.31; 0 -0.06 0 0.06; -0.06 0 0.06 0; 1 1 1 1; 0 0 0 0; 0.37 0.37 0.37 0.37; 0 0 0 0; 0.06 0.06 0.06 0.06]; 
%Geometria que gera a superficie minimima que é um circulo de raio 0.06 e centrado em (0.37, 0)
geometryFromEdges(model,dl);
generateMesh(model,'Hmax',0.006);
%pdemesh(model); axis equal;

% Estimar os termos da qe J_\phi - ordens
% Checar/colocar a mash na posição original
% Campo mgnetico gerado pela corrente


%%
aux1 = model.Mesh.Nodes(1,:); %meus xizes
aux2 = model.Mesh.Nodes(2,:); %meus yipsilons
CAMPO=load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65_2.mat');
G.out=campo2(aux1,aux2); %sai os campos em função dos xizes e yipsilons
if 0 figure(8); pdemesh(model); axis equal; end
G.out1 = resistivity_nova(G.n0,G.Te0,1);
v_ei = G.out1.v_ei;  G.rho = G.out1.eta_par; G.constante = 7.89e11*(1.38e-23)*300;
v_en= G.constante*(G.ng-G.n0); G.somatranpost = v_ei + v_en;

applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',@inifuncadapt);
setInitialConditions(model,@inifuncadapt); 

specifyCoefficients(model,'m',0,'d',@d,'c',@c,'a',@a,'f',@s03_Gaussiana);

%resolver o sistema
results = solvepde(model,tlist);
%plotar a solução
sol=results.NodalSolution;
%plotasol2(model,sol)
%figure(7); pdeplot(model,'XYData',G.plott1);  colormap(jet); title('v_i_o_n-v_l_o_s_s'); axis equal
if 0
%gera_Anim(model,sol,[1 1])
pause(1)
gera_Anim(model,sol,[1 2])
pause(1)
%gera_Anim(model,sol,[1 3])
%pause(1)
%gera_Anim(model,sol,[1 4])
%pause(1)
end
%plotasol_anim(model,sol)
C=[1 2 3 4];
%plotasolTermo(model,sol)
%
plotasol(model,sol,C)
keyboard
end

