function Simulacao_n_J()
clear all
close all
clc
global G e me mi CAMPO
Neq=4; e= 1.6e-19; me = 9.11e-31; mi = 1.67e-27;
Ntime = 1000;            indice=8; % numero que sera salvo os resultados da simulação
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
G.ligas=0; % Ativa a fonte apoximada pela gaussiana
%%
% u(1,:) -> n
% u(2,:) -> J_phi
% u(3,:) -> pe
% u(4,:) -> pion
model = createpde(Neq);  % dl = decsg(gd,sf,ns)  pdecirc(0.37,0,0.06)  pdecirc(0.37,0,0.01)
dl=[1 1 1 1; 0.31 0.37 0.43 0.37; 0.37 0.43 0.37 0.31; 0 -0.06 0 0.06; -0.06 0 0.06 0; 1 1 1 1; 0 0 0 0; 0.37 0.37 0.37 0.37; 0 0 0 0; 0.06 0.06 0.06 0.06]; 
%Geometria que gera a superficie minimima que é um circulo de raio 0.06 e centrado em (0.37, 0)
geometryFromEdges(model,dl);
generateMesh(model,'Hmax',0.003);
%pdemesh(model); axis equal;

% Estimar os termos da qe J_\phi - ordens
% Checar/colocar a mash na posição original
% Campo mgnetico gerado pela corrente


%%
aux1 = model.Mesh.Nodes(1,:); %meus xizes
aux2 = model.Mesh.Nodes(2,:); %meus yipsilons
CAMPO=load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65_2.mat');
G.out=campo2(aux1,aux2); %sai os campos em função dos xizes e yipsilons
if 0
figure(5); pdeplot(model,'XYData',sqrt(G.out.Br.^2 +G.out.Bz.^2)); title('Campo magnético externo poloidal'); axis(G.eixo); colormap('jet'); axis equal;
figure(6); pdeplot(model,'XYData',G.out.Ephi); title('Campo elétrico externo toroidal'); axis(G.eixo); colormap('jet'); axis equal;
figure(7); pdeplot(model,'XYData',G.out.Bphi); title('Campo magnético externo toroidal'); axis(G.eixo); colormap('jet'); axis equal;
figure(8); pdemesh(model); axis equal;
%if G.ligacorrente figure(9); pdeplot(model,'XYData',sqrt(G.out.Bply.^2 +G.out.Bplx.^2)); title('Campo magnético poloidal plasma corrente unitaria'); axis(G.eixo); colormap('jet'); axis equal; end
end

G.out1 = resistivity_nova(G.n0,G.Te0,1);
v_ei = G.out1.v_ei;  G.rho = G.out1.eta_par; G.constante = 7.89e11*(1.38e-23)*300;
v_en= G.constante*(G.ng-G.n0); G.somatranpost = v_ei + v_en;


applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',@inifuncadapt01);
setInitialConditions(model,@inifuncadapt01); 
%@c01
specifyCoefficients(model,'m',0,'d',@d01,'c',@c01,'a',@a01,'f',@s03);

%resolver o sistema
results = solvepde(model,tlist);
%plotar a solução
sol=results.NodalSolution;
%plotasol2(model,sol)
%keyboard
G.plotton=1;
%figure(20); pdeplot(model,'XYData',G.plott10);  colormap(jet); title('v_i_o_n-v_l_o_s_s'); axis equal
if 0
G.plotton=1;
%gera_Anim(model,sol,[1 1]); pause(1)
gera_Anim(model,sol,[1 2]) ;pause(1) %gera_Anim(model,sol,[1 3]); pause(1); gera_Anim(model,sol,[1 4]); pause(1);
end
%keyboard
if 0 %anima em prop o J_phi  
    figure(12); s=length(sol(3,2,:)); p=round(s*0.15/15,0); T=1:p:(s-1); 
    amax=max(max(sol(:,1,:)./sol(:,1,1)));
        aux1=0.026;%sol(:,3,1)./sol(:,1,1)/e;
    amax2=max(max(sol(:,3,:)./sol(:,1,:)/e/aux1));
    for oll=T
        aux=sol(:,3,oll)./sol(:,1,oll)/e;
        subplot(1,2,1); pdeplot(model,'XYData',aux/aux1); X2=['T_e(:,:,dt*' num2str(oll) ')/T_e_0']; box on; title(X2); colormap(jet); axis equal; caxis([1 amax2]);
        subplot(1,2,2); pdeplot(model,'XYData',sol(:,1,oll)./sol(:,1,1)); X2=['n(:,:,dt*' num2str(oll) ')/n_0']; box on; title(X2); colormap(jet); axis equal; caxis([1 amax]); pause(0.10);
    end
end
C=[1 2 3 4];
%plotasolTermo(model,sol)
%plotasol(model,sol,C)
%plotasol_prop(model,sol,C)
%keyboard
if 1
          %plot proporcional
         figure(1);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',G.plottv0-G.plottl0); 
         box on;  colormap(jet);  axis equal; 
         X=['(v_i_n-v_l_o_s_s)(:,1)']; title(X);
         
         subplot(1,2,2);  pdeplot(model,'XYData',G.plottv6-G.plottl6);
         box on;  colormap(jet); axis equal; 
         X=['(v_i_n-v_l_o_s_s)(:,60)']; title(X); 
         ti2=['vinvloss1']; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico

         
         it = [Ntime/8 2*Ntime/8 3*Ntime/8 4*Ntime/8 5*Ntime/8 6*Ntime/8 7*Ntime/8 Ntime];
         %plot da densidade numerica
         i=1; amax=max(max(sol(:,i,:)./sol(:,i,1))); 
         figure(1);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(1))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['n(:,:,' num2str(it(1)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(2))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['n(:,:,' num2str(it(2)) ')']; title(ti1); 
         ti2=['ntod1B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         %%%%%%%%%%
         figure(2);   set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(3))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['n(:,:,' num2str(it(3)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(4))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['n(:,:,' num2str(it(4)) ')']; title(ti1); 
         ti2=['ntod2B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         %%%%%%%%%%
         figure(3);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(5))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['n(:,:,' num2str(it(5)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(6))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['n(:,:,' num2str(it(6)) ')']; title(ti1); 
         ti2=['ntod3B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         
         figure(4);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(7))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['n(:,:,' num2str(it(7)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(8))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['n(:,:,' num2str(it(8)) ')']; title(ti1); 
         ti2=['ntod4B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         
               %plot da densidade de corrente
         i=2; amax=max(max(sol(:,i,:)./sol(:,i,1))); 
         figure(1);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(1))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['Jphi(:,:,' num2str(it(1)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(2))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['Jphi(:,:,' num2str(it(2)) ')']; title(ti1); 
         ti2=['Jphitod1B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         %%%%%%%%%%
         figure(2);   set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(3))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['Jphi(:,:,' num2str(it(3)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(4))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['Jphi(:,:,' num2str(it(4)) ')']; title(ti1); 
         ti2=['Jphitod2B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         %%%%%%%%%%
         figure(3);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(5))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['Jphi(:,:,' num2str(it(5)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(6))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['Jphi(:,:,' num2str(it(6)) ')']; title(ti1); 
         ti2=['Jphitod3B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
                  %%%%%%%%%%
         figure(4);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(7))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['Jphi(:,:,' num2str(it(7)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(8))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['Jphi(:,:,' num2str(it(8)) ')']; title(ti1); 
         ti2=['Jphitod4B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         
               %plot da pressao eletronica
         i=3; amax=max(max(sol(:,i,:)./sol(:,i,1))); 
         figure(1);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(1))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pe(:,:,' num2str(it(1)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(2))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pe(:,:,' num2str(it(2)) ')']; title(ti1); 
         ti2=['petod1B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         %%%%%%%%%%
         figure(2);   set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(3))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pe(:,:,' num2str(it(3)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(4))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pe(:,:,' num2str(it(4)) ')']; title(ti1); 
         ti2=['petod2B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         %%%%%%%%%%
         figure(3);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(5))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pe(:,:,' num2str(it(5)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(6))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pe(:,:,' num2str(it(6)) ')']; title(ti1); 
         ti2=['petod3B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         %%%%%%%%%%
         figure(4);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(7))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pe(:,:,' num2str(it(7)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(8))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pe(:,:,' num2str(it(8)) ')']; title(ti1); 
         ti2=['petod4B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         
           
               %plot da pressao ionica 
         i=4; amax=max(max(sol(:,i,:)./sol(:,i,1))); 
         figure(1);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(1))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pi(:,:,' num2str(it(1)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(2))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pi(:,:,' num2str(it(2)) ')']; title(ti1); 
         ti2=['pitod1B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         %%%%%%%%%%
         figure(2);   set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(3))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pi(:,:,' num2str(it(3)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(4))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pi(:,:,' num2str(it(4)) ')']; title(ti1); 
         ti2=['pitod2B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
         %%%%%%%%%%
         figure(3);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(5))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pi(:,:,' num2str(it(5)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(6))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pi(:,:,' num2str(it(6)) ')']; title(ti1); 
         ti2=['pitod3B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
          %%%%%%%%%%
         figure(4);  set(gcf, 'Position',  [1, 401, 1200, 400]); 
         subplot(1,2,1);  pdeplot(model,'XYData',sol(:,i,it(7))./sol(:,i,1)); 
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pi(:,:,' num2str(it(7)) ')']; title(ti1); 
         
         subplot(1,2,2);  pdeplot(model,'XYData',sol(:,i,it(8))./sol(:,i,1));
         box on;  colormap(jet); caxis([1 amax]); axis equal; 
         ti1=['pi(:,:,' num2str(it(8)) ')']; title(ti1); 
         ti2=['pitod4B' num2str(indice)]; ti = [ti2 '.png']; saveas(gcf,ti);          %salvar grafico
                 
    
end
%keyboarde
end
