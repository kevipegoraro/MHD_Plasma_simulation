c = 1;
a = 0;
f = @circlef;
numberOfPDE = 1;
model = createpde(numberOfPDE);
g = @circleg;
geometryFromEdges(model,g);
figure; 
pdegplot(model,'EdgeLabels','on'); 
axis equal
title 'Geometry With Edge Labels Displayed';
applyBoundaryCondition(model,'dirichlet','Edge',(1:4),'u',0);
[u,p,e,t] = adaptmesh(g,model,c,a,f,'tripick','circlepick','maxt',2000,'par',1e-3);
figure; 
pdemesh(p,e,t); 
axis equal
x = p(1,:)';
y = p(2,:)';
r = sqrt(x.^2+y.^2);
uu = -log(r)/2/pi;
figure;
pdeplot(p,e,t,'XYData',u-uu,'ZData',u-uu,'Mesh','off');
figure;
pdeplot(p,e,t,'XYData',u,'ZData',u,'Mesh','off');