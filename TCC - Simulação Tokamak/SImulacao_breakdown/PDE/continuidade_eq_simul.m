%clear all; close all; clc
global G e me mi
Neq=1; e= 1.6e-19; me = 9.11e-31; mi = 1.67e-27;
Ntime = 3000;            indice=8; % numero que sera salvo os resultados da simulação
dt = 1e-5; tlist = 0:dt:dt*Ntime;
G.R0 = 0.37;  G.control=0;  G.t=0;
G.B0 = 0.7*2;   G.p0 = 0.05; 
G.ng = G.p0/1.38e-23/300;  G.n0 = 1e8;
G.Te0 = 0.026;             G.gas = 'H2'; 
G.Vloop = 10;              G.a0=0.06;  G.eixo=[G.R0-G.a0 G.R0+G.a0 -G.a0 G.a0];
G.C1 = 1;                   G.V = 1;     G.camp=1;  
G.ligacorrente = 1; % Ativa o calculo do campo magnetico gerado pelo plasma
G.ligacoef=1;
G.cplott = 0; G.plotton = 1; % Ativa o plot da fonte em todos os tempos
G.ligas=0.001; % Ativa a fonte apoximada pela gaussiana
%%
% u(1,:) -> n
model = createpde(Neq);  % dl = decsg(gd,sf,ns)  pdecirc(0.37,0,0.06)  pdecirc(0.37,0,0.01)
dl=[1 1 1 1; 0.31 0.37 0.43 0.37; 0.37 0.43 0.37 0.31; 0 -0.06 0 0.06; -0.06 0 0.06 0; 1 1 1 1; 0 0 0 0; 0.37 0.37 0.37 0.37; 0 0 0 0; 0.06 0.06 0.06 0.06]; 
%Geometria que gera a superficie minimima que é um circulo de raio 0.06 e centrado em (0.37, 0)
geometryFromEdges(model,dl);
generateMesh(model,'Hmax',0.003);
%pdemesh(model); axis equal;
%%
G.out1 = resistivity_nova(G.n0,G.Te0,1);
v_ei = G.out1.v_ei;  G.rho = G.out1.eta_par; G.constante = 7.89e11*(1.38e-23)*300;
v_en= G.constante*(G.ng-G.n0); G.somatranpost = v_ei + v_en;
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',G.n0);
setInitialConditions(model,G.n0); 
D=0.04; %state.time*
f = @(location,state)state.time*10e10*exp(-1000*((-G.R0+location.x).^2+(location.y+0.04-state.time*3).^2))+...
    (0.02-state.time)*10e10*exp(-1000*((-G.R0+location.x).^2+(location.y-0.04+state.time*3).^2))+...
    (0.02-state.time)*10e10*exp(-1000*((-G.R0+location.x+0.04-state.time*3).^2+(location.y).^2))+...
    state.time*10e10*exp(-1000*((-G.R0+location.x-0.04+state.time*3).^2+(location.y).^2));
c = @(location,state)D*location.x;
d = @(location,state)location.x;
specifyCoefficients(model,'m',0,'d',d,'c',c,'a',0,'f',f);

%resolver o sistema
results = solvepde(model,tlist);
%plotar a solução
sol=results.NodalSolution;
%keyboard
s=length(sol(1,:));
if 1 %anima em prop o J_phi  
    figure(2); p=15;%round(s*0.15/15,0); 
    T=1:p:Ntime; 
    amax=max(max(sol(:,:)./sol(:,1)));
    for oll=T
        aux=sol(:,oll)./sol(:,1);
        pdeplot(model,'XYData',aux); X=['n(:,']; X2=[X 'dt*' num2str(oll) ')']; box on; title(X2); colormap(jet); axis equal; caxis([1 amax]);
        pause(0.07);
    end
end
keyboard
if 0
 X=['n(:,'];
T=[3 round(s/4,0) round(3*s/4,0) s-1];
amax=max(max(sol(:,:)./sol(:,1)));
subplot(2,2,1)
pdeplot(model,'XYData',sol(:,T(1))./sol(:,1))
X2=[X 'dt*' num2str(T(1)) ')'];
title(X2)
colormap(jet); caxis([1 amax]);
axis equal 
subplot(2,2,2)
pdeplot(model,'XYData',sol(:,T(2))./sol(:,1))
X2=[X 'dt*' num2str(T(2)) ')'];
title(X2)
colormap(jet); caxis([1 amax]);
axis equal 
subplot(2,2,3)
pdeplot(model,'XYData',sol(:,T(3))./sol(:,1))
X2=[X 'dt*' num2str(T(3)) ')'];
title(X2)
colormap(jet); caxis([1 amax]);
axis equal 
subplot(2,2,4)
pdeplot(model,'XYData',sol(:,T(4))./sol(:,1))
X2=[X 'dt*' num2str(T(4)) ')'];
title(X2)
colormap(jet); caxis([1 amax]);
axis equal
end
