
%%
close all
clear all
D=0.001;
model = createpde(1);
g = @circleg;
geometryFromEdges(model,g);
applyBoundaryCondition(model,'dirichlet','Edge',(1:4),'u',0);
u0 = @u0fun;
setInitialConditions(model,u0);
ff=@s01;
specifyCoefficients(model,'m',0,'d',1,'c',1,'a',0,'f',1);
keyboard
%%
generateMesh(model,'Hmax',0.05);
tlist = 0:0.01:0.5;
results = solvepde(model,tlist);
sol = results.NodalSolution;

subplot(2,2,1)
pdeplot(model,'XYData',sol(:,3))
title('Time 0.02')
axis equal 
subplot(2,2,2)
pdeplot(model,'XYData',sol(:,5))
title('Time 0.04')
axis equal 
subplot(2,2,3)
pdeplot(model,'XYData',sol(:,11))
title('Time 0.1')
axis equal 
subplot(2,2,4)
pdeplot(model,'XYData',sol(:,51))
title('Time 0.5')
axis equal 
%%

setInitialConditions(model,results)
tlist1 = 0.5:0.01:1.0;
results1 = solvepde(model,tlist1);
%%
sol1 = results1.NodalSolution;
figure
subplot(2,2,1)
pdeplot(model,'XYData',sol1(:,1))
title('Time 0.5')
axis equal 
subplot(2,2,2)
pdeplot(model,'XYData',sol1(:,21))
title('Time 0.7')
axis equal 
subplot(2,2,3)
pdeplot(model,'XYData',sol1(:,41))
title('Time 0.9')
axis equal 
subplot(2,2,4)
pdeplot(model,'XYData',sol1(:,51))
title('Time 1.0')
axis equal 
%%

setInitialConditions(model,results,21)
%%

tlist2 = 0.2:0.01:1.0;
results2 = solvepde(model,tlist2);

