%Solve system 2 equations
model = createpde(2);
geometryFromEdges(model,@circleg);
specifyCoefficients(model,'m',0,...
                          'd',[1;1],...
                          'c',[1;0],...
                          'a',[0;1],...
                          'f',@s02)
%%
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',[0 0]);
setInitialConditions(model,0);
generateMesh(model,'Hmax',0.25);
tlist = 0:0.1:5;
results = solvepde(model,tlist);
sol=results.NodalSolution;
%%
% View the solution.
%pdeplot(model,'XYData',sol)
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
%keyboard
axis equal
keyboard

%lala