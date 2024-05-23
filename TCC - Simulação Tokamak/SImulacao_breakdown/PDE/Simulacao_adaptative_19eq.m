%%
global out e me mi Neq nrpontos
Neq=10; e= 1.6e-19; me = 9.11e-31; mi = 1.67e-27;
Ntime = 5; dt = 1e-2; tlist = 0:dt:dt*Ntime;
R0 = 0.37; B0 = 0.7;  Vloop = 10;  
p0 = 0.05; %0.01 % 1e-4 milibar = 0.01 Pa
ng = p0/1.38e-23/300; %particulas do gás
n0 = 1e8; %densidade de particulas inicial do plasma
Te0 = 0.026; % teperatura inicial do gas
gas = 'H2'; %tipo de gas
D=10;

%%
% u(1,:) -> n
% u(2,:) -> J_x
% u(3,:) -> J_y
% u(4,:) -> J_phi
% u(5,:) -> pe
% u(6,:) -> pion
% u(7,:) -> v_en
% u(8,:) -> v_in
% u(9,:) -> v_ei
% u(10,:) -> v_loss
% u(11,:) -> B_pl_x
% u(12,:) -> B_pl_y
% u(13,:) -> B_pl_phi
% u(14,:) -> E_pl_x
% u(15,:) -> E_pl_y
% u(16,:) -> E_pl_phi
% u(17,:) -> nabla X E_pl_x 
% u(18,:) -> nabla X E_pl_x
% u(19,:) -> nabla X E_pl_x

f = @circlef;
model = createpde(Neq); pini = createpde(1); 
g = @circleg;
geometryFromEdges(model,g);
%figure; pdegplot(model,'EdgeLabels','on'); axis equal; title 'Geometry With Edge Labels Displayed';
applyBoundaryCondition(model,'dirichlet','Edge',(1:4),'u',0);
[u,p,ei,t] = adaptmesh(g,pini,1,0,f,'tripick','circlepick','maxt',1000,'par',1e-3);
%% Plot Finest Mesh
%figure; pdemesh(p,ei,t); axis equal
%keyboard
setInitialConditions(model,@inifunc); 
%condição constante na fronteira do circulo
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0); 
%Cada coluna p (:, k) consiste na coordenada x do ponto k em p (1, k) e na coordenada y do ponto k em p (2, k).
out=campo(R0+p(1,:)/10,p(2,:)/10);
nrpontos = size(out.Br,2);
%Plot do campo
%Bpol = sqrt(out.Br.^2 +out.Bz.^2); pdeplot(p,ei,t,'XYData',Bpol)

%%
specifyCoefficients(model,'m',0,...
                          'd',d,... 
                          'c',c,...
                          'a',a,...
                          'f',@s03_adapt)

%%

results = solvepde(model,tlist);
%plotar a solução
sol=results.NodalSolution;
figure;
pdeplot(p,ei,t,'XYData',u,'ZData',u,'Mesh','off');

