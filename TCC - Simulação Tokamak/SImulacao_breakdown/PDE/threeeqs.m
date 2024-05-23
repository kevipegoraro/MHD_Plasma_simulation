model = createpde(2);
geometryFromEdges(model,@circleg);
u0 = @u0fun;
ut0 = @ut0fun;
generateMesh(model,'Hmax',0.25);
setInitialConditions(model,0,0);
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',[1 0]);
%setInitialConditions(model,10); %[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
specifyCoefficients(model,'m',0,...
                          'd',[1 ; 1],...
                          'c',[10; 1],...
                          'a',0,...
                          'f',@s02)
tlist = 0:0.01:0.5;
results = solvepde(model,tlist);
%setInitialConditions(model,results)
%results = solvepde(model,tlist);
sol=results.NodalSolution;
size(sol)
%%
% View the solution.
%pdeplot(model,'XYData',sol)
figure(1)
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

figure(2)
subplot(2,2,1)
pdeplot(model,'XYData',sol(:,2,3))
title('Time 0.02')
axis equal 
subplot(2,2,2)
pdeplot(model,'XYData',sol(:,2,5))
title('Time 0.04')
axis equal 
subplot(2,2,3)
pdeplot(model,'XYData',sol(:,2,11))
title('Time 0.1')
axis equal 
subplot(2,2,4)
pdeplot(model,'XYData',sol(:,2,51))
title('Time 0.5')
%keyboard
axis equal
keyboard