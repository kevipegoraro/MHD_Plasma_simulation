% Rotina principal do PDE solver
global out e me mi Neq rho out1 ng a1 n0 Te0 Leff R0 B0 p0 constante v_ei  somatranpost  Vloop % Declaracao de variaveis globais que serao usadas em outros scripts
Neq=4; % Numero de equacoes do sistema 
e= 1.6e-19; me = 9.11e-31; mi = 1.67e-27; % Constantes fisicas 
R0 = 0.37; B0 = 0.7;  Vloop = 10; p0 = 0.05; % Dados do Tokamak NOVA-FURG  
Ntime = 20;  dt = 1e-4; tlist = 0:dt:dt*Ntime; % Tempos onde a solucao sera calculada
ng = p0/1.38e-23/300; %particulas do gas
n0 = 1e8; %densidade de particulas inicial do plasma
Te0 = 0.026; % teperatura inicial do gas
gas = 'H2'; %tipo de gas
%%
% Quem sera cada elemento de u
% u(1,:) -> n
% u(2,:) -> J_phi
% u(3,:) -> pe
% u(4,:) -> pion

model = createpde(Neq);  % Criando o sistema de EDP's
geometryFromEdges(model,@circleg); % Criando uma geometria ciscular
generateMesh(model,'Hmax',0.15); % Gerando a mesh com precisão 0.15
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',@inifuncadapt); % Aplicando as condicoes de borda
setInitialConditions(model,@inifuncadapt);  % Aplicando as condicoes iniciais

 %%
aux1 = model.Mesh.Nodes(1,:); % Vetor de localizacoes na malha na direcao x
aux2 = model.Mesh.Nodes(2,:); % Vetor de localizacoes na malha na direcao y
out=campo2(aux1,aux2); % Retorna os campos eletricos e magneticos externos em função de aux1 e aux2

pdemesh(model); axis equal; title('Mesh triangular usada'); axis([-1 1 -1 1]) % Plot da mesh
pdeplot(model,'XYData',sqrt(out.Br.^2 +out.Bz.^2)); title('Campo magnetico externo poloidal'); axis([-1 1 -1 1]); colormap('jet'); % Plot do campo magnetico poloidal

out1 = resistivity_nova(n0,Te0,1); % Retorna uma estrutura com o coeficiente de tranpsote de ionizacao eletron-ion e eta.
v_ei = out1.v_ei; % Definindo o coeficiente de tranpsote de ionizacao eletron-ion
a1 = first_townsend_coeff(p0,out.Ephi,gas,0); % Retorna o primeiro doeficiente de townsend
Leff = R0*B0./(0.001+(out.Br.^2 + out.Bz.^2).^0.5); % Retorna a distancia media percorida pelas particulas antes de colidir com a parede
prop=0.6; % Fator de convercao do raio da secao de vacuo do tokamak para o raio unitario da simulacao
rho = out1.eta_par; % Retorna o eta
constante = 7.89e11*(1.38e-23)*300; v_en=constante*(ng-n0); % Definicao do coeficiente de ionizacao eletron neutros
somatranpost = v_ei + v_en; % Somando o coeficiente de ionizacao eletron-ion e o coeficiente de ionizacao eletron neutros
D=2; % Definindo o coeficiente de difusao
d=ones(Neq,1); % Definindo o coeficiente d das derivadas lineares de cada equacao do sistema
c=zeros(Neq,1); c(1,1)=D;  % Definindo o coeficiente do termo do laplaciano, que esta apenas na equacao da continuidade
specifyCoefficients(model,'m',0,'d',d,'c',c,'a',@a,'f',@s03_adapt_simpli) % Fixando os coeficientes
% Resolver o sistema
results = solvepde(model,tlist);
% Plotar a solucao
sol=results.NodalSolution;
plotasol2(model,sol)
%fim

