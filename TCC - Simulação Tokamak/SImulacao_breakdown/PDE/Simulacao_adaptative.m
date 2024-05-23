clear all
%close all
clc
global out e me mi Neq rho out1 ng a1 n0 Te0 Leff R0 B0 p0 constante v_ei 
Neq=4; e= 1.6e-19; me = 9.11e-31; mi = 1.67e-27;
Ntime = 20; 
dt = 1e-4; tlist = 0:dt:dt*Ntime;
R0 = 0.37; 
B0 = 0.7;  Vloop = 10;  
p0 = 0.05; %0.01 % 1e-4 milibar = 0.01 Pa
ng = p0/1.38e-23/300; %particulas do gás
n0 = 1e8; %densidade de particulas inicial do plasma
Te0 = 0.026; % teperatura inicial do gas
gas = 'H2'; %tipo de gas

%%
% u(1,:) -> n
% u(2,:) -> J_phi
% u(3,:) -> pe
% u(4,:) -> pion
model = createpde(Neq); 
geometryFromEdges(model,@circleg);
generateMesh(model,'Hmax',0.2);
%%
%[u,p,ei,t] = adaptmesh(g,model,1,0,f,'tripick','circlepick','maxt',50,'par',1e-3);
%% Plot Finest Mesh
%figure; pdemesh(p,ei,t); axis equal
%keyboard
%%
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',@inifuncadapt);
setInitialConditions(model,@inifuncadapt); 
%%
%condição constante na fronteira do circulo
%applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0); 
%%
aux1 = model.Mesh.Nodes(1,:); %x
aux2 = model.Mesh.Nodes(2,:); %y
aux=2*length(aux1);
A=zeros(2,aux);
A(1,1) = aux1(1); A(1,2) = aux1(1);
A(2,1) = aux2(1); A(2,2) = aux2(1);
for j=3:1:aux
    ai=round(j/2,0);
    A(1,j)= (aux1(ai)+aux1(ai-1))/2;
    A(2,j)= (aux2(ai)+aux2(ai-1))/2;
end

out=campo(R0+A(1,:)*0.06,A(2,:)*0.06);
%Plot do campo Bpol = sqrt(out.Br.^2 +out.Bz.^2); pdeplot(model,'XYData',Bpol)
out1 = resistivity_nova(n0,Te0,1);
v_ei = out1.v_ei;
a1 = first_townsend_coeff(p0,out.Ephi,gas,0);
Leff = R0*B0./(0.001+(out.Br.^2 + out.Bz.^2).^0.5);
rho = out1.eta_par;
constante = 7.89e11*(1.38e-23)*300;
v_en=constante*(ng-n0);
somatranpost = v_ei + v_en;
cPe = -2*n0*(e^2)*rho/mi;
cPi = cPe;
%%
D=0.15;
d=ones(Neq,1); c=zeros(Neq,1); a=c; 
a(2,1) = somatranpost; a(3,1) = cPe; a(4,1) = cPi;
c(1,1)=D;
specifyCoefficients(model,'m',0,'d',d,'c',c,'a',a,'f',@s03_adapt)
%resolver o sistema
results = solvepde(model,tlist);
%plotar a solução
sol=results.NodalSolution;
%pdeplot(p,ei,t,'XYData',sol(:,1,2))%'ZData',u,'Mesh','off');
C=[1 2 3];
plotasol(model,sol,C)
%keyboard
