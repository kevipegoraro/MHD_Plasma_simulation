function Simulacao_TCABR()
clear all
close all
clc
global G e me mi 
Neq=4; e= 1.6e-19; me = 9.11e-31; mi = 1.67e-27;
Ntime = 1000;            indice=8; % numero que sera salvo os resultados da simulação
dt = 1e-5; tlist = 0:dt:dt*Ntime;
G.R0 = 0.37;  G.control=0;  G.t=0;
G.B0 = 0.7*2;   G.p0 = 0.05; 
G.ng = G.p0/1.38e-23/300;  G.n0 = 1e8;
G.Te0 = 0.026;             G.gas = 'H2'; 
G.Vloop = 10;              G.a0=0.06;  G.eixo=[G.R0-G.a0 G.R0+G.a0 -G.a0 G.a0];
G.C1 = 1; G.V = 1; G.camp=1; G.ligacorrente = 0; G.ligacoef=0; G.cplott = 0; G.plotton = 1; 
G.ligas=1; 
%%
% u(1,:) -> n
% u(2,:) -> J_phi
% u(3,:) -> pe
% u(4,:) -> pion
%keyboard
model = createpde(Neq);  % dl = decsg(gd,sf,ns)  pdecirc(0.37,0,0.06)  pdecirc(0.37,0,0.01)
dl= [
    2.0000    2.0000    2.0000    2.0000
    0.8500    0.3950    0.3950    0.8500
    0.3950    0.3950    0.8500    0.8500
   -0.2720   -0.2720    0.2720    0.2720
   -0.2720    0.2720    0.2720   -0.2720
         0         0         0         0
    1.0000    1.0000    1.0000    1.0000] ;
geometryFromEdges(model,dl);
generateMesh(model,'Hmax',0.01);
%pdemesh(model); axis equal;
G.out1 = resistivity_nova(G.n0,G.Te0,1);
v_ei = G.out1.v_ei;  G.rho = G.out1.eta_par; G.constante = 7.89e11*(1.38e-23)*300;
v_en= G.constante*(G.ng-G.n0); G.somatranpost = v_ei + v_en;
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',@ini_TCABR);
setInitialConditions(model,@ini_TCABR); 
%@c01 @a01_teste [-0.02;0;0;0]  [0;0.01;0;0]
specifyCoefficients(model,'m',0,'d',@d_TCABR,'c',@c_TCABR,'a',@a_TCABR,'f',@f_TCABR);

results = solvepde(model,tlist);
sol=results.NodalSolution;
%plotasol_prop(model,sol,[1 2 3 4])
if 1 
     p=25;%round(s*0.15/15,0); 
    T=1:p:Ntime; 
    amax=max(max(sol(:,1,:)./sol(:,1,1)));
    amax2=max(max(-sol(:,2,:)./sol(:,2,1)));
    figure(1); set(gcf, 'Position',  [1, 800, 1000, 450]);
    %figure(2); set(gcf, 'Position',  [600, 600, 600, 450]);
    %keyboard
    for oll=T
        aux=sol(:,1,oll)./sol(:,1,1); %figure(1); 
        subplot(1,2,1); pdeplot(model,'XYData',aux); X=['n(:,']; X2=[X 'dt*' num2str(oll) ')']; box on; title(X2); colormap(jet); axis equal; caxis([1 amax]);
        pause(0.3);  %figure(2);
        subplot(1,2,2); pdeplot(model,'XYData',sol(:,2,oll)./sol(:,2,1)); X=['j(:,']; X2=[X 'dt*' num2str(oll) ')']; box on; title(X2); colormap(jet); axis equal; 
        %pause(0.3);
    end
end

end
