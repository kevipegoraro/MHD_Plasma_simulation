% - teste malha adaptativa
%%
Ntime = 5; dt = 1e-5; R0 = 0.37; n0 = 1e2; Te0 = 0.026; gas = 'H2'; a0=0.07;
%n=n0*ones(tam,Ntime+1); %definindo a condição inicial da densidade de particulas
x0=a0/4; sgma=-10; s0=1; D=0.001;
c = 1;
a = 0;
f = @circlef;
numberOfPDE = 1;
model = createpde(numberOfPDE);
g = @circleg;
geometryFromEdges(model,g);

%figure; 
%pdegplot(model,'EdgeLabels','on'); 
%axis equal


applyBoundaryCondition(model,'dirichlet','Edge',(1:4),'u',0);
[u,p,e,t] = adaptmesh(g,model,c,a,f,'tripick','circlepick','maxt',2000,'par',1e-3);
%%
tamanho=size(u)
x = p(1,:)';
y = p(2,:)';
r = sqrt((x+0.3).^2+y.^2);
size(r)
s=s0*exp(r/sgma);
nr = length(region.x); 
f = zeros(nr); 
f = region.x - region.y + state.u(1,:);

specifyCoefficients(model,'m',0,...
                          'd',1,...
                          'c',D,...
                          'a',0,...
                          'f',f);

%%
generateMesh(model,'Hmax',0.05);
tlist = 0:0.01:0.5;
results = solvepde(model,tlist);
% Plot the solution for times 0.02, 0.04, 0.1, and 0.5.
sol = results.NodalSolution;

subplot(2,2,1)
pdeplot(model,'XYData',sol(:,3))
title('Time 0.02')
subplot(2,2,2)
pdeplot(model,'XYData',sol(:,5))
title('Time 0.04')
subplot(2,2,3)
pdeplot(model,'XYData',sol(:,11))
title('Time 0.1')
subplot(2,2,4)
pdeplot(model,'XYData',sol(:,51))
title('Time 0.5')

% Plot the 2-D solution.
% pdeplot(model,'XYData',u)
% % figure(1);
% % pdeplot(p,e,t,'XYData',uu,'ZData',uu,'Mesh','off');
% % figure(2);
% % pdeplot(p,e,t,'XYData',uu);    
% % %pdeplot(p,e,t,'XYData',uu);
% axis equal
% 

%%

% figure; %Plote os valores de erro.
% pdeplot(p,e,t,'XYData',u-uu,'ZData',u-uu,'Mesh','off');
% figure;
% pdeplot(p,e,t,'XYData',u,'ZData',u,'Mesh','off');
%keyboard